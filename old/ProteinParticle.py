#%% Import all packages in need.
from enum import Enum
import random
import math
import Parameters as pa
#%% Define the Enum class for different types of proteins.
class Protein(Enum):
    FREE=0
    BOUNDARY=1
    DEFAULT=2
    LIGHT_INDUCED=3
    DEFAULT_FLUORESCENT=-2
    LIGHT_INDUCED_FLUORESCENT=-3
    def get_fluorescent_type(self):
        if (self in [Protein.FREE,Protein.BOUNDARY]):
            return self
        else:
            return Protein(-abs(self.value))
    def get_nonfluorescent_type(self):
        if (self in [Protein.FREE,Protein.BOUNDARY]):
            return self
        else:
            return Protein(abs(self.value))
    def is_fluorescent(self):
        '''
        Determine whether the particle is fluorescent.
        If it is free or boundary, return -1.
        If it is fluorescent, return 1.
        If it is non-fluorescent, return 0.
        '''
        if (self in [Protein.FREE,Protein.BOUNDARY]):
            return -1
        elif (self.value>0):
            return 0
        else:
            return 1

if (len(list(vars(Protein).items()))-11!=pa.numProteinType):#Check whether the number or types of protein is consistent with parameters.
    raise ValueError('''
          The number of protein types is inconsistent with the parameters!
    ''')
#%% Define the operations for the numbers and vectors indicating different directions.    
dictDirection={0:(0,0,0),1:(-1,0,0),2:(1,0,0),3:(0,0,1),4:(0,0,-1),5:(0,1,0),6:(0,-1,0)}#The correspondence between direction numbering and the coordinate change.
dictOppositeDirection={0:0,1:2,2:1,3:4,4:3,5:6,6:5}#The correspondence between a direction and its opposite direction.
def plus_tuple(x,y:tuple):
    x_=list(x)
    x_[0]+=y[0]
    x_[1]+=y[1]
    x_[2]+=y[2]
    if type(x)==tuple:
        return tuple(x_)
    else:
        return x_
#%% Define the class for the grid.
class Grid:
    grid={}
    bx=0
    by=0
    bz=0
    def __init__(self,bx:int,by:int,bz:int):
        '''
        The function to initialize a grid with the given size.
        bx, by and bz represent the size of the particles' range of motion.
        All positions of the grid, except its boundary, are filled with "Protein.FREE".
        '''
        print('''
              Initializing the grid...
        ''')
        self.bx=bx
        self.by=by
        self.bz=bz
        for i in range(1,bx):
            for j in range(1,by):
                for k in range(1,bz):
                    self.grid[(i,j,k)]=Protein.FREE
        for i in [0,bx]:
            for j in range(by+1):
                for k in range(bz+1):
                    self.grid[(i,j,k)]=Protein.BOUNDARY
        for j in [0,by]:
            for i in range(bx+1):
                for k in range(bz+1):
                    self.grid[(i,j,k)]=Protein.BOUNDARY
        for k in [0,bz]:
            for i in range(bx+1):
                for j in range(by+1):
                    self.grid[(i,j,k)]=Protein.BOUNDARY
        return None
    def num_surrounding(self,x:int,y:int,z:int,proteinType:Protein):
        '''
        This function counts the number of protein with the type "proteinType" within the surrounding 26 positions of (x,y,z).
        '''
        num=0
        if (x<=0 or x>=self.bx or y<=0 or y>=self.by or z<=0 or z>=self.bz):
            num=26
        else:
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    for k in [-1,0,1]:
                        if my_grid.grid[plus_tuple((x,y,z),(i,j,k))]==proteinType:
                            num+=1
            if my_grid.grid[(x,y,z)]==proteinType:
                num-=1
        return num

my_grid=Grid(*pa.b)
listProteinParticles=[]
#%% Define the class for single protein particles.
class Particle:
    '''
    x: The particle's corrdinate in the three-dimensional grid.
    v: The velocity of the particle.
    t: The types of protein particles within the surrounding of this particle, in the order of itself, left, right, up, down, forward, back.
    p: The probability that the particle stays, moves to the left, right, up, down, forward, back.
    '''
    x=[0,0,0]
    v=[0,0,0]
    t=[Protein.DEFAULT,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE]
    p=[1,0,0,0,0,0,0]
    def __init__(self,x:int,y:int,z:int,vx:float,vy:float,vz:float,t:'Protein'=Protein.DEFAULT):
        if (x<=0 or x>=my_grid.bx or y<=0 or y>=my_grid.by or z<=0 or z>=my_grid.bz):
            raise ValueError
        self.x=[x,y,z]
        self.v=[vx,vy,vz]
        self.t=[t,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE,Protein.FREE]
        return None
    def velocity_square(self):
        return self.v[0]*self.v[0]+self.v[1]*self.v[1]+self.v[2]*self.v[2]
    def swap(self,particle2:'Particle'):
        '''
        Swap the properties of two particles.
        '''
        tmp=Particle(*tuple(self.x),*tuple(self.v),self.t[0])
        tmp.t=self.t
        self.x=particle2.x
        self.v=particle2.v
        self.t[0]=particle2.t[0]
        particle2.x=tmp.x
        particle2.v=tmp.v
        particle2.t[0]=tmp.t[0]
        my_grid.grid[tuple(self.x)]=self.t[0]
        my_grid.grid[tuple(particle2.x)]=particle2.t[0]
    def renew_t(self,grid:dict):
        '''
        Renew the surrounding protein types stored in "t" based on the information of the grid.
        '''
        self.t[1:]=[grid[plus_tuple(tuple(self.x),dictDirection[i])] for i in range(1,7)]
    def allocate_p(self):
        '''
        Allocate the probabilities that the protein particle moves in all directions based on its velocity and the types of surrounding protein particles.
        '''
        self.renew_t(my_grid.grid)
        if self.v[0]>0:
            self.p[1]=0
            self.p[2]=self.v[0]*self.v[0]/self.velocity_square()*(1-math.exp(-pa.Lambda*self.v[0]))
        else:
            self.p[1]=self.v[0]*self.v[0]/self.velocity_square()*(1-math.exp(pa.Lambda*self.v[0]))
            self.p[2]=0
        if self.v[2]>0:
            self.p[3]=self.v[2]*self.v[2]/self.velocity_square()*(1-math.exp(-pa.Lambda*self.v[2]))
            self.p[4]=0
        else:
            self.p[3]=0
            self.p[4]=self.v[2]*self.v[2]/self.velocity_square()*(1-math.exp(pa.Lambda*self.v[2]))
        if self.v[1]>0:
            self.p[5]=self.v[1]*self.v[1]/self.velocity_square()*(1-math.exp(-pa.Lambda*self.v[1]))
            self.p[6]=0
        else:
            self.p[5]=0
            self.p[6]=self.v[1]*self.v[1]/self.velocity_square()*(1-math.exp(pa.Lambda*self.v[1]))
        self.p[0]=1-sum(self.p[1:])
        for direction in range(1,7):
            if self.t[direction] not in [Protein.FREE,Protein.BOUNDARY]:
                self.p[direction]*=pa.e
        for proteinNameIndexThis in range(pa.numProteinType):
            if Protein[pa.listProteinNames[proteinNameIndexThis]]==self.t[0]:
                break
        for direction in range(1,7):
            for proteinNameIndex in range(2,pa.numProteinType):
                proteinType=Protein[pa.listProteinNames[proteinNameIndex]]
                numSurroundingNow=my_grid.num_surrounding(*tuple(self.x),proteinType)
                numSurroundingNext=my_grid.num_surrounding(*plus_tuple(tuple(self.x),dictDirection[direction]),proteinType)
                if numSurroundingNext<numSurroundingNow:
                    reduceFactor=math.exp((numSurroundingNext-numSurroundingNow)*pa.Delta_E[proteinNameIndex][proteinNameIndexThis]/pa.k_B/pa.T)
                    self.p[0]+=self.p[direction]*(1-reduceFactor)
                    self.p[direction]*=reduceFactor
        return None
    def reallocate_p(self):
        '''
        Reallocate the probabilities that the protein particle moves in all directions such that the probabilities adds up to 1.
        The ratio between the probability the protein stays and moves should be equal to the ratio of the two parameters of the function.
        The relative proportion of the probability of moving in each direction does not change.
        '''
        self.allocate_p()
        p_stay=self.p[0]
        p_move=sum(self.p[1:])
        if abs(p_stay<1e-5):
            p_stay=0
        if abs(p_move<1e-5):
            p_move=0
        if p_stay<0 or p_move<0:
            raise ValueError('''
                 'p_stay' and 'p_move' must be non-negative!
                             ''')
        if p_stay==0 and p_move==0:
            self.p=[1/7]*7
        else:
            self.p[0]=p_stay/(p_stay+p_move)
            sum_p_move=sum(self.p[1:])
            if sum_p_move==0:
                for i in range(1,7):
                    self.p[i]=p_move/(p_stay+p_move)/6
            else:
                for i in range(1,7):
                    self.p[i]=p_move/(p_stay+p_move)*self.p[i]/sum_p_move
        return None
    def direction_determine(self):
        '''
        Generate a random number to determine which direction to move in.
        Return the random number generated an integer as a parameter for the "move" function.
        '''
        self.reallocate_p()
        r=random.uniform(0,1)
        for direction in range(6):
            if r<=sum(self.p[:direction+1]):
                return r,direction
        return r,6
    def move(self,direction:int):
        '''
        The value of the parameter "direction" represents:
        0: stay
        1: left
        2: right
        3: up
        4: down
        5: forward
        6: back
        '''
        if direction not in [0,1,2,3,4,5,6]:
            raise ValueError
        if direction!=0:
            previous_x=self.x
            next_x=plus_tuple(self.x,dictDirection[direction])
            if my_grid.grid[tuple(next_x)]==Protein.BOUNDARY:
                #next_x=[(next_x[i]-1)%(pa.b[i]-1)+1 for i in range(3)]#Cyclic boundary condition.
                next_x=plus_tuple(previous_x,dictDirection[dictOppositeDirection[direction]])#Reflective boundary condition.
            if my_grid.grid[tuple(next_x)]==Protein.FREE:
                self.x=[next_x[i] for i in range(3)]
                for i in range(3):
                    if self.v[i]*dictDirection[direction][i]>0:
                        self.v[i]*=pa.f
                my_grid.grid[tuple(previous_x)]=Protein.FREE
                my_grid.grid[tuple(self.x)]=self.t[0]
            else:
                for particle in listProteinParticles:
                    if particle.x==next_x:
                        self.swap(particle)
                        break
        return self
    def fluorescent(self,f:bool):
        '''
        Set the fluorescent state of this particle.
        If "f" is "True", set the particle to fluorescent.
        If "f" is "False", set the particle to fluorescent.
        '''
        if f:
            self.t[0]=self.t[0].get_fluorescent_type()
        else:
            self.t[0]=self.t[0].get_nonfluorescent_type()
        my_grid.grid[tuple(self.x)]=self.t[0]
#%% Define the function to do fluorescence quenching
def quench(listProteinParticles:list):
    pass
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes

def str_func(funcString:str,ls:list):
    if funcString.lower()=='max':
        return max(ls)
    elif funcString.lower()=='min':
        return min(ls)
    elif funcString.lower()=='sum':
        return sum(ls)
    elif funcString.lower()=='avg':
        return sum(ls)/len(ls)
    else:
        raise ValueError('Invalid function!')

def projection(listProteinParticles:list,my_grid:'Grid',unit_size:tuple,proj_direction:str,proj_type:str):
    '''
    This function takes the average fluorescent intensity of each region of given size in the grid, and project it to a 2D plane.
    listProteinParticles: the list of all protein particles.
    grid: the dictionary of the type of protein particles in each grid.
    unit_size: the size of each region of which to take the average fluorescent intensity. Should be a tuple with two entries.
    proj_direction: the direction of projecting. Should be one of 'x', 'y' or 'z'.
    proj_type: decide how to show the projection, to use the maximal or average value of all planes. Should be one of 'avg', 'max', 'min' or 'sum'. Each protein particle is assigned a value, 1 for fluorescent and 0 for non-fluorescent.
    Return a 2D array, each entry for the value of its corresponding unit.
    '''
    if len(unit_size)!=2 or abs(int(unit_size[0]-1))!=unit_size[0]-1 or abs(int(unit_size[1]-1))!=unit_size[1]-1:
        raise ValueError(''''unit_size' should be a two-entry tuple of positive integers.''')
    direction_index={'x':0,'y':1,'z':2}[proj_direction]
    grid_size=(my_grid.bx-1,my_grid.by-1,my_grid.bz-1)
    unit_number=(int(grid_size[direction_index-2]/unit_size[0]),int(grid_size[direction_index-1]/unit_size[1]))
    output=[]
    for i in range(unit_number[0]):
        output.append([])
        for j in range(unit_number[1]):
            sum_plane=[0]*grid_size[direction_index]
            for k in range(1,grid_size[direction_index]):
                for l in range(i*unit_size[0]+1,(i+1)*unit_size[0]+1):
                    for m in range(j*unit_size[1]+1,(j+1)*unit_size[1]+1):
                        position=[0,0,0]
                        position[direction_index]=k
                        position[direction_index-2]=l
                        position[direction_index-1]=m
                        position=tuple(position)
                        if my_grid.grid[position].is_fluorescent()==1:
                            sum_plane[k]+=1
            output[i].append(str_func(proj_type,sum_plane))
    return output

def visualize(output:list,vmin=0,vmax=0,saveFig:bool=False,savePath:str=''):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    if vmin>=vmax:
        im=ax.imshow(output,cmap='BuGn')
    else:
        im=ax.imshow(output,cmap='BuGn',vmin=vmin,vmax=vmax)
    plt.colorbar(im)
    plt.axis('off')
    if saveFig:
        if savePath=='':
            raise AttributeError('''
                 When 'saveFig' is True, 'savePath' cannot be empty!
            ''')
        plt.savefig(savePath)
    else:
        plt.show()
    plt.close()
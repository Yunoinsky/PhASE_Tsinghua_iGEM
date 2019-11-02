#ifndef _PROTEIN_PARTICLE_H
#define _PROTEIN_PARTICLE_H
#include <math.h>
#include <vector>
#include "parameters.h"
#include "frequently_used_functions.h"
class Grid{
    const int bx,by,bz;
    class Particle ****grid;
public:
    Grid(int,int,int);
    ~Grid();
    int get_size(int);
    enum Protein position(std::vector<int>);
    enum Protein position(int,int,int);
    void point_to(std::vector<int>,const class Particle&);
    void point_to(int,int,int,const class Particle&);
    void point_to(std::vector<int>);
    void point_to(int,int,int);
    class Particle* get_ref(std::vector<int>);
    class Particle* get_ref(int,int,int);
    bool in_interior(int,int,int);
}my_grid(pa::b[0],pa::b[1],pa::b[2]);
bool Grid::in_interior(int x,int y,int z){
    if(x<=0||x>=bx||y<=0||y>=by||z<=0||z>=bz) return false;
    else return true;
}
Grid::Grid(int bx,int by,int bz):bx(bx),by(by),bz(bz){
    /*
    The function to initialize a grid with the given size.
    bx, by and bz represent the size of the particles' range of motion.
    All positions of the grid, except its boundary, are filled with "FREE".
    */
    printf("Initializing the grid...\n");
    grid=(Particle****)malloc((bx+1)*sizeof(Particle***));
    for(int i=0;i<bx+1;i++){
        grid[i]=(Particle***)malloc((by+1)*sizeof(Particle**));
        for(int j=0;j<by+1;j++){
            grid[i][j]=(Particle**)malloc((bz+1)*sizeof(Particle*));
            for(int k=0;k<bz+1;k++) grid[i][j][k]=NULL;
        }
    }
    printf("The grid has been initialized!\n");
}
Grid::~Grid(){
    for(int i=bx;i>=0;i--){
        for(int j=by;j>=0;j--) free(grid[i][j]);
        free(grid[i]);
    }
    free(grid);
}
int Grid::get_size(int index){
    int sz;
    switch(index){
        case 0:sz=bx;break;
        case 1:sz=by;break;
        case 2:sz=bz;break;
    }
    return sz;
}
void Grid::point_to(std::vector<int> index,const Particle& p){
    point_to(index[0],index[1],index[2],p);
}
void Grid::point_to(int x,int y,int z,const Particle& p){
    grid[x][y][z]=const_cast<Particle *>(&p);
}
void Grid::point_to(std::vector<int> index){
    point_to(index[0],index[1],index[2]);
}
void Grid::point_to(int x,int y,int z){
    grid[x][y][z]=NULL;
}
Particle* Grid::get_ref(std::vector<int> index){
    return get_ref(index[0],index[1],index[2]);
}
Particle* Grid::get_ref(int x,int y,int z){
    if(x<0||x>bx||y<0||y>by||z<0||z>bz) return NULL;
    return grid[x][y][z];
};
const int numParticle=sum(std::vector<int>(pa::N+2,pa::N+pa::numProteinType));

class Particle{
    std::vector<int> x;//The particle's corrdinate in the three-dimensional grid.
    std::vector<double> v;//The velocity of the particle.
    std::vector<enum Protein> t;//The types of protein particles within the surrounding of this particle, in the order of itself, left, right, up, down, forward, back.
    std::vector<double> p;//The probability that the particle stays, moves to the left, right, up, down, forward, back.
public:
    Particle(int,int,int,double,double,double,enum Protein);
    int get_x(int index){
        return x[index];
    }
    enum Protein get_t(int index){
        return t[index];
    }
    double get_v(int index){
        return v[index];
    }
    double velocity_square() const{
        double vsq=0;
        for(int i=0;i<3;i++) vsq+=v[i]*v[i];
        return vsq;
    }
    void add_v(int index,double v_incre){
        v[index]+=v_incre;
    }
    void swap(Particle&);
    void renew_t();
    void allocate_p();
    enum Direction direction_determine();
    void move(enum Direction);
    void set_fluorescent(bool);
};
int num_surrounding(int x,int y,int z,enum Protein proteinType){
    int num=0;
    if(!(my_grid.in_interior(x,y,z))) num=0;
    else{
        for(int i=-1;i<=1;i++)
            for(int j=-1;j<=1;j++)
                for(int k=-1;k<=1;k++)
                    if(my_grid.get_ref(x+i,y+j,z+k)!=NULL&&my_grid.get_ref(x+i,y+j,z+k)->get_t(0)==proteinType) num++;      
        if(my_grid.get_ref(x,y,z)!=NULL&&my_grid.get_ref(x,y,z)->get_t(0)==proteinType) num--;
    }
    return num;
}
int num_surrounding(std::vector<int> index,enum Protein proteinType){
    return num_surrounding(index[0],index[1],index[2],proteinType);
}
Particle::Particle(int cx,int cy,int cz,double vx,double vy,double vz,enum Protein pt):x(3),v(3),t(7),p(7){
    if(!my_grid.in_interior(cx,cy,cz)&&pt!=BOUNDARY){
        printf("The coordinate (%d, %d, %d) is out of range!\n",cx,cy,cz);
        exit(1);
    }
    else if(pt==BOUNDARY&&!my_grid.in_interior(cx,cy,cz)){
        if(cx<0) cx=0;
        else if(cx>my_grid.get_size(0)) cx=my_grid.get_size(0);
        if(cy<0) cy=0;
        else if(cy>my_grid.get_size(1)) cy=my_grid.get_size(1);
        if(cz<0) cz=0;
        else if(cz>my_grid.get_size(2)) cz=my_grid.get_size(2);        
    }
    int x[3]={cx,cy,cz};
    double v[3]={vx,vy,vz};
    enum Protein t[7]={pt,FREE,FREE,FREE,FREE,FREE,FREE};
    double p[7]={1,0,0,0,0,0,0};
    this->x=std::vector<int>(x,x+3);
    this->v=std::vector<double>(v,v+3);
    this->t=std::vector<enum Protein>(t,t+7);
    this->p=std::vector<double>(p,p+7);
    my_grid.point_to(this->x,*this);
}
void Particle::swap(Particle& p2){
    /*
    Swap the properties of two particles.
    */
    Particle tmp(x[0],x[1],x[2],v[0],v[1],v[2],t[0]);
    tmp.t=t;
    x=p2.x;
    v=p2.v;
    t[0]=p2.t[0];
    p2.x=tmp.x;
    p2.v=tmp.v;
    p2.t[0]=tmp.t[0];
    my_grid.point_to(x,*this);
    my_grid.point_to(p2.x,p2);
}
void Particle::renew_t(){
    /*
    Renew the surrounding protein types stored in "t" based on the information of the grid.
    */
    for(int i=1;i<7;i++) t[i]=my_grid.position(x);
}
void Particle::allocate_p(){
    /*
    Allocate the probabilities that the protein particle moves in all directions based on its velocity and the types of surrounding protein particles.
    */
    renew_t();
    int p_index[6]={1,2,6,5,4,3};
    for(int i=0;i<3;i++){
        if(v[i]>0){
            p[p_index[2*i]]=0;
            p[p_index[2*i+1]]=v[i]*v[i]/velocity_square()*(1-exp(-pa::Lambda*v[i]));
        }
        else{
            p[p_index[2*i]]=v[i]*v[i]/velocity_square()*(1-exp(pa::Lambda*v[i]));
            p[p_index[2*i+1]]=0;
        }
    }
    p[0]=1-sum(p,1,7);
    for(int i=1;i<7;i++){
        if(t[i]!=FREE&&t[i]!=BOUNDARY) p[i]*=pa::e;
    }
    double reduceFactor;
    // double reduceFactorModifier;
    // Particle *reduceFactorNeighbor;
    int indexProteinTypeThis=int(t[0]);
    for(int i=1;i<7;i++){
        for(int j=2;j<pa::numProteinType;j++){
            enum Protein proteinType=Protein(j);
            if(pa::N[j]==0) continue;
            int numSurroundingNow=num_surrounding(x,proteinType);
            std::vector<int> moveDirectionVector(direction_vector(Direction(i)));
            int numSurroundingNext=num_surrounding(x+moveDirectionVector,proteinType);
            if(numSurroundingNext<numSurroundingNow){
                reduceFactor=exp((numSurroundingNext-numSurroundingNow)*(*(pa::Delta_E+j*pa::numProteinType+indexProteinTypeThis))/pa::R/pa::T);
                // reduceFactorNeighbor=my_grid.get_ref(x+direction_vector(opposite(Direction(i))));
                // if(reduceFactorNeighbor!=NULL){
                //     if(reduceFactorNeighbor->p[i]>0) reduceFactorModifier=reduceFactorNeighbor->p[i];
                //     else reduceFactorModifier=0;
                // }
                // else reduceFactorModifier=0;
                // reduceFactor=reduceFactor*(1-reduceFactorModifier)+reduceFactorModifier*pa::e;
                p[0]+=p[i]*(1-reduceFactor);
                p[i]*=reduceFactor;
            }
        }
    }
    double p_stay=p[0],p_move=sum(p,1,7);
    if(abs(p_stay)<1e-5) p_stay=0;
    if(abs(p_move)<1e-5) p_move=0;
    if(p_stay<0||p_move<0){
        printf("p_stay and p_move must be non-negative!\n");
        exit(1);
    }
    double p_total=p_stay+p_move;
    if(p_total==0)
        for(int i=0;i<7;i++) p[i]=0.1428571429;
    else
        for(int i=0;i<7;i++) p[i]/=p_total;
}
enum Direction Particle::direction_determine(){
    /*
    Generate a random number to determine which direction to move in.
        Return the random number generated an integer as a parameter for the "move" function.
    */
    allocate_p();
    double r=random_uniform(0,1);
    for(int i=0;i<6;i++){
        if(r<=p[i]) return Direction(i);
        else r-=p[i];
    }
    return Direction(6);
}
void Particle::move(enum Direction direction){
    if(direction!=ORIGIN){
        std::vector<int> previous_x(x),next_x(x+direction_vector(direction));
        if(my_grid.position(next_x)==BOUNDARY){
            //for(int i=0;i<3;i++) next_x[i]=(next_x[i]+pa::b[i]-2)%(pa::b[i]-1)+1;//Cyclic boundary condition.
            next_x=previous_x+direction_vector(opposite(direction));//Reflective boundary condition.
        }
        else if(my_grid.position(next_x)==FREE){
            x=next_x;
            for(int i=0;i<3;i++)
                if(v[i]*direction_vector(direction)[i]>0) v[i]*=pa::f;
            my_grid.point_to(previous_x);
            my_grid.point_to(x,*this);
        }
        else swap(*my_grid.get_ref(next_x));
    }
}
void Particle::set_fluorescent(bool f){
    /*
    Set the fluorescent state of this particle.
    If "f" is "True", set the particle to fluorescent.
    If "f" is "False", set the particle to fluorescent.
    */
    if(f) t[0]=fluorescent_type(t[0]);
    else t[0]=non_fluorescent_type(t[0]);
}
enum Protein Grid::position(std::vector<int> index){
    return position(index[0],index[1],index[2]);
}
enum Protein Grid::position(int x,int y,int z){
    if(x<0||x>bx||y<0||y>by||z<0||z>bz) return FREE;
    else if(grid[x][y][z]==NULL){
        if(in_interior(x,y,z)) return FREE;
        else return BOUNDARY;
    }
    else return grid[x][y][z]->get_t(0);
}
#endif
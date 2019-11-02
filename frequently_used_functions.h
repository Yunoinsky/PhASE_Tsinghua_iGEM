#ifndef _FREQUENTLY_USED_FUNCTIONS_H
#define _FREQUENTLY_USED_FUNCTIONS_H
#include <stdio.h>
#ifdef WIN32
    #include <winbase.h>
#else
    #include <sys/time.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <mm_malloc.h>
#include <math.h>
#include <vector>
const double PI=3.141592654;

/*Protein enum type.*/
    enum Direction{ORIGIN=0,LEFT=1,RIGHT=2,UP=3,DOWN=4,FORWARD=5,BACK=6};
    enum Direction opposite(enum Direction d){
        enum Direction o;
        switch(d){
            case ORIGIN:o=ORIGIN;break;
            case LEFT:o=RIGHT;break;
            case RIGHT:o=LEFT;break;
            case UP:o=DOWN;break;
            case DOWN:o=UP;break;
            case FORWARD:o=BACK;break;
            case BACK:o=FORWARD;break;
        }
        return o;
    }
    std::vector<int> direction_vector(enum Direction d){
        std::vector<int> dv(3);
        for(int i=0;i<3;i++) dv[i]=0;
        switch(d){
            case ORIGIN:break;
            case LEFT:dv[0]=-1;break;
            case RIGHT:dv[0]=1;break;
            case UP:dv[2]=1;break;
            case DOWN:dv[2]=-1;break;
            case FORWARD:dv[1]=1;break;
            case BACK:dv[1]=-1;break;
        }
        return dv;
    }

    enum Protein{FREE=0,BOUNDARY=1,DEFAULT=2,LIGHT_INDUCED=3,DEFAULT_FLUORESCENT=4,LIGHT_INDUCED_FLUORESCENT=5};
    enum Protein fluorescent_type(enum Protein n){
        enum Protein f;
        switch(n){
            case FREE:
            case BOUNDARY:
            case DEFAULT_FLUORESCENT:
            case LIGHT_INDUCED_FLUORESCENT:f=n;break;
            case DEFAULT:f=DEFAULT_FLUORESCENT;break;
            case LIGHT_INDUCED:f=LIGHT_INDUCED_FLUORESCENT;break;
        }
        return f;
    }
    enum Protein non_fluorescent_type(enum Protein f){
        enum Protein n;
        switch(f){
            case FREE:
            case BOUNDARY:
            case DEFAULT:
            case LIGHT_INDUCED:n=f;break;
            case DEFAULT_FLUORESCENT:n=DEFAULT;break;
            case LIGHT_INDUCED_FLUORESCENT:n=LIGHT_INDUCED;break;
        }
        return n;
    }
    bool is_fluorescent(enum Protein p){
        /*
        Determine whether the particle is fluorescent.
        */
        bool i;
        switch(p){
            case FREE:
            case BOUNDARY:
            case DEFAULT:
            case LIGHT_INDUCED:i=false;break;
            case DEFAULT_FLUORESCENT:
            case LIGHT_INDUCED_FLUORESCENT:i=true;break;
        }
        return i;
    }
/*End*/

/*Tuple (array) operations.*/
    template <typename T>
    std::vector<T> operator+ (std::vector<T> a,std::vector<T> b){
        size_t n=a.size();
        if(n!=b.size()){
            printf("Two vectors to add are not of the same size!\n");
            exit(1);
        }
        else{
            std::vector<T> result(n);
            for(size_t i=0;i<n;i++) result[i]=a[i]+b[i];
            return result;
        }
    }
    template <typename T>
    bool operator== (std::vector<T> a,std::vector<T> b){
        bool result=true;
        size_t n=a.size();
        if(n!=b.size()) return false;
        for(size_t i=0;i<n;i++)
            if(a[i]!=b[i]){
                result=false;
                break;
            }
        return result;
    }
    template <typename T>
    T sum(std::vector<T> a,size_t start_index,size_t end_index){
        T result=0;
        for(size_t i=start_index;i<end_index;i++) result+=a[i];
        return result;
    }
    template <typename T>
    T sum(std::vector<T> a){
        size_t start_index=0,end_index=a.size();
        return sum(a,start_index,end_index);
    }
    std::vector<int> num_to_position(size_t num,int bx,int by,int bz){
        std::vector<int> position(3);
        if(num<=0||num>(bx-1)*(by-1)*(bz-1)){
            printf("Number out of range in function \"num_to_position\"");
            exit(1);
        }
        size_t tmp_num=num;
        position[0]=(tmp_num-1)%(bx-1)+1;
        tmp_num=(tmp_num-position[0])/(bx-1);
        position[1]=tmp_num%(by-1)+1;
        tmp_num=(tmp_num-position[1]+1)/(by-1);
        position[2]=tmp_num+1;
        return position;
    }
    size_t position_to_num(std::vector<int> position,int bx,int by,int bz){
        return (size_t)((bx-1)*(by-1)*(position[2]-1)+(bx-1)*(position[1]-1)+position[0]);
    }
/*End*/

/*Random number generation.*/
    double random_uniform(double a,double b,long int *seed){
        double t;
        *seed=2045*(*seed)+1;
        *seed=*seed-(*seed/1048576)*1048576;
        t=(*seed)/1048576.0;
        t=a+(b-a)*t;
        return t;
    }
    double random_uniform(double a,double b){
        srand(clock());
        long int s=rand();
        return random_uniform(a,b,&s);
    }
    double random_gauss(double  mu, double sigma){
        static double U,V;
        static int phase=0;
        double z;
        srand(clock());
        if(phase == 0){
            U=rand()/(RAND_MAX+1.0);
            V=rand()/(RAND_MAX+1.0);
            z=sqrt(-2.0*log(U))*sin(2.0*PI*V);
        }
        else{
            z=sqrt(-2.0*log(U))*cos(2.0*PI*V);
        }
        phase=1-phase;
        return z*sigma+mu;
    }
    std::vector<int> n_choose_m(int n,int m) {
        std::vector<int> result(n);
        result.clear();
        result.reserve(n);
        srand(clock());
        for(size_t i=0;i<n;i++){
            result.push_back(i);
        }
        int p1;
        int p2;
        int temp;
        while (--n){
            p1 = n;
            p2 = rand() % n;
            temp = result[p1];
            result[p1] = result[p2];
            result[p2] = temp;
        }
        std::vector<int> result_m(m);
        for(size_t i=0;i<m;i++) result_m[i]=result[i]+1;
        return result_m;
    }
/*End*/

/*File Operation*/
    void copy_file(const char* ori,const char* dest){
        FILE *op, *inp;
        op=fopen(ori,"r");
        char filePath[36];
        sprintf(filePath,"%s%s%s",dest,"/",ori);
        inp=fopen(filePath,"w");
        void *buf;
        while(!feof(op)){
            fread(&buf,1,1,op);
            fwrite(&buf,1,1,inp);
        }
        fclose(inp);
        fclose(op);
    }
/*End*/
#endif
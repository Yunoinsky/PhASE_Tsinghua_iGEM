#include <stdio.h>
#include <string.h>
#include <mm_malloc.h>
#include <math.h>
#include <limits.h>
#include <vector>
#include "parameters.h"
#include "frequently_used_functions.h"
#include "protein_particle.h"
#ifdef WIN32
    void create_folder(const char* folderName){
        system(strcat("md ",folderName));
    }
#else
    #include <sys/stat.h>
    #include <sys/types.h>
    void create_folder(const char* folderName){
        mkdir(folderName,S_IRWXU|S_IRWXG|S_IRWXO);
    }
#endif

const bool inputFromFile=false;//If input from file, the file name should be "command.txt".
const int input_batch=5;//The maximum number of characters to read from the input.
int string_to_input(const char* c){
    int result;
    switch(c[0]){
        case 'f':
        case 'F':result=-1;break;
        case 'p':
        case 'P':result=-2;break;
        case 'j':
        case 'J':result=-3;break;
        case 'l':
        case 'L':result=-4;break;
        case 'w':
        case 'W':result=-5;break;
        case 'v':
        case 'V':result=-6;break;
        case 'r':
        case 'R':result=-7;break;
        case 'q':
        case 'Q':result=INT_MIN;break;
        default:{
            result=0;
            for(int i=0;i<input_batch&&c[i]>='0'&&c[i]<='9';i++) result=result*10+int(c[i]-'0');
        }
    }
    return result;
}
int input_from_file(){
    char inputLine[input_batch];
    static FILE *fptr=NULL;
    if(fptr==NULL&&(fptr=fopen("command.txt","r"))==NULL){
        printf("Run into error when opening the command file!\n");
        exit(1);
    }
    fgets(inputLine,4,fptr);
    return string_to_input(inputLine);
}
int input_from_console(){
    printf("Type a positive integer (less than 100000) as the iteration time.\nType \"f\" to start the FRAP experiment.\nType \"p\" to print the current image of the system.\nType \"j\" to print the projection image of the system.\nType \"l\" to add light induction.\nType \"w\" to withdraw light induction.\nType \"v\" to create video files.\nType \"r\" to recover from log file.\nType \"q\" to quit the program.\n");
    char inputLine[input_batch];
    scanf("%s",inputLine);
    return string_to_input(inputLine);
}
int (*get_input)()=inputFromFile?input_from_file:input_from_console;

int main(int argc,char* argv[]){
    std::vector<int> listPositionNumbers=n_choose_m((pa::b[0]-1)*(pa::b[1]-1)*(pa::b[2]-1),numParticle);
    std::vector<std::vector<int> > listPositions(numParticle,std::vector<int>(3));
    std::vector<Particle> listParticles(0,Particle(0,0,0,0,0,0,BOUNDARY));
    int indexProteinType=2;
    int numTypeProtein=0;
    double random_v[3],std_v;
    for(int i=0;i<numParticle;i++){
        listPositions[i]=num_to_position(listPositionNumbers[i],pa::b[0],pa::b[1],pa::b[2]);
        numTypeProtein++;
        if(numTypeProtein>pa::N[indexProteinType]){
            numTypeProtein=0;
            indexProteinType++;
            i--;
            continue;
        }
        std_v=sqrt(pa::k_B*pa::T/pa::m[indexProteinType]);
        for(int j=0;j<3;j++) random_v[j]=random_gauss(0,std_v);
        listParticles.push_back(Particle(listPositions[i][0],listPositions[i][1],listPositions[i][2],random_v[0],random_v[1],random_v[2],Protein(indexProteinType)));
    }
    char folderPath[23];
    time_t timep;
    time(&timep);
    strftime(folderPath,sizeof(folderPath),"OutputFigures-%m%d%H%M",localtime(&timep));
    create_folder(folderPath);
    int inputNumber=0,countRun=0;
    char fileName[14],filePath[23];
    FILE *logFile;
    copy_file("parameters.h",folderPath);
    while(inputNumber!=INT_MIN){
        inputNumber=get_input();
        while(inputNumber>0){
            countRun++;
            sprintf(fileName,"/%07d.csv",countRun);
            strcpy(filePath,folderPath);
            strcat(filePath,fileName);
            logFile=fopen(filePath,"w");
            fprintf(logFile,"%s","Particle Address,Type,x[0],x[1],x[2],v[0],v[1],v[2],Direction Determined\n");
            for(int i=0;i<numParticle;i++){
                fprintf(logFile,"%p,%d,%d,%d,%d,%e,%e,%e",&listParticles[i],int(listParticles[i].get_t(0)),listParticles[i].get_x(0),listParticles[i].get_x(1),listParticles[i].get_x(2),listParticles[i].get_v(0),listParticles[i].get_v(1),listParticles[i].get_v(2));
                int directionDetermined=listParticles[i].direction_determine();
                fprintf(logFile,",%d\n",directionDetermined);
                listParticles[i].move(Direction(directionDetermined));
            }
            std::vector<double> v_square_avg(pa::numProteinType,0);
            for(int i=0;i<numParticle;i++){
                v_square_avg[int(listParticles[i].get_t(0))]+=listParticles[i].velocity_square();
            }
            for(int numTypeProtein=2;numTypeProtein<pa::numProteinType;numTypeProtein++){
                if(pa::N[numTypeProtein]>0) {
                    v_square_avg[numTypeProtein]/=pa::N[numTypeProtein];
                    if(v_square_avg[numTypeProtein]<pa::k_B*pa::T/pa::m[indexProteinType]*3){
                        double std_v_incr=sqrt(pa::k_B*pa::T/pa::m[indexProteinType]-v_square_avg[numTypeProtein]/3);
                        for(int i=0;i<numParticle;i++){
                            if(int(listParticles[i].get_t(0))==numTypeProtein){
                                for (int j=0;j<3;j++) listParticles[i].add_v(j,random_gauss(0,std_v_incr));
                            }
                        }
                    }
                }

            }
            fclose(logFile);
            inputNumber--;
        }
        if(inputNumber==0) continue;
        else if(inputNumber==-4) pa::Delta_E=&pa::Delta_E_light[0][0];
        else if(inputNumber==-5) pa::Delta_E=&pa::Delta_E_nolight[0][0];
    }
    printf("The program exits normally.\n");
    return 0;
}
#ifndef _PARAMETERS_H
#define _PARAMETERS_H
namespace pa{
    const int numProteinType=6;
//Hyperparameters
    const double e=2;//The factor to increase the probability of protein-protein swapping.
    const double Lambda=1e7;
    double Delta_E_nolight[numProteinType][numProteinType]={
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,1e-2,1e-2,1e-2,1e-2},
        {0,0,1e-2,1e-2,1e-2,1e-2},
        {0,0,1e-2,1e-2,1e-2,1e-2},
        {0,0,1e-2,1e-2,1e-2,1e-2}
    };//The (i,j)-th entry represents the binding energy between the i-th protein and the j-th protein.
    double Delta_E_light[numProteinType][numProteinType]={
        {0,0,0,0,0,0},
        {0,0,0,0,0,0},
        {0,0,1e-2,1e-2,1e-2,1e-2},
        {0,0,1e-2,1e-2,1e-2,1e-2},
        {0,0,1e-2,1e-2,1e10,1e10},
        {0,0,1e-2,1e-2,1e10,1e10}
    };//The (i,j)-th entry represents the binding energy between the i-th protein and the j-th protein.
    const double *Delta_E=&Delta_E_light[0][0];//The (i,j)-th entry represents the binding energy between the i-th protein and the j-th protein.
    const double m[numProteinType]={0,0,6e-10,6e-10,6e-10,6e-10};//Masses of all types of protein particles.
    const double f=0.5;//The factor to reduce the velocity at the direction in which the protein moves because of the resistance.
//Experimental conditions
    const int b[3]={200,100,100};//The size of the particles' range of motion.
    const double T=298;//Temperature.
    const int N[numProteinType]={0,0,0,0,0,5000};//The initial number of each kind of protein particle.
//Scientific constants
    const double R=8.314;//Ideal gas constant.
    const double k_B=1.380649e-23;//Boltzmann constant.
}
#endif
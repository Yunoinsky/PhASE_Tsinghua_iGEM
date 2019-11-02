#%% Import parameters from the Excel workbook.
import numpy as np
import pandas as pd
from openpyxl import load_workbook
parameterFile='ProteinProperties.xlsx'
df=pd.read_excel(parameterFile,sheet_name='Input')
numProteinType=len(df)-3
listProteinNames=list(df.PPI[0:numProteinType])
#%% Hyperparameters
e=2#The factor to increase the probability of protein-protein swapping.
Lambda=1e7
Delta_E=[list(map(float,row)) for row in np.array(df.iloc[0:numProteinType,1:numProteinType+1]).tolist()]#The (i,j)-th entry represents the binding energy between the i-th protein and the j-th protein.
m=list(map(float,np.array(df.iloc[numProteinType+2,1:numProteinType+1]).tolist()))#Masses of all types of protein particles. 
f=0.5#The factor to reduce the velocity at the direction in which the protein moves because of the resistance.
#%% Experimental conditions
b=(200,100,100)#The size of the particles' range of motion.
T=298#Temperature.
N=list(map(float,np.array(df.iloc[numProteinType+1,1:numProteinType+1]).tolist()))#The initial number of each kind of protein particle.
projUnitSize=(2,2)#The unit size of projection.
#%% Scientific constants
R=8.314#Ideal gas constant.
k_B=1.380649e-23#Boltzmann constant.
#%% A function to record parameters.
def record_parameters(excelFileName:str):
    book=load_workbook(excelFileName)
    writer=pd.ExcelWriter(excelFileName,engine ='openpyxl') 
    writer.book=book 
    parameterIndex=['Hyper parameters','e','Lambda','f','Experimental conditions','bx','by','bz','T/K','Projection unit size','Scientific constants','R/(J`mol^-1`K^-1)','k_B/(J`K^-1)']
    parameterValue=['',e,Lambda,f,'',b[0],b[1],b[2],T,projUnitSize,'',R,k_B]
    df=pd.DataFrame(parameterValue,index=parameterIndex,columns=['Value'])
    df.to_excel(writer,sheet_name='Parameters')
    writer.save()
    writer.close()
#%% Run self-test.
raiseValueError=False
raiseValueErrorMessage=''
if (e<1):
    raiseValueErrorMessage+='''
          "e" must be greater than or equal to 1!
    '''
    raiseValueError=True
if (Lambda<=0):
    raiseValueErrorMessage+='''
          "Lambda" must be greater than 0!
    '''
    raiseValueError=True
numProteinType=len(Delta_E)
for row in Delta_E:
    isSymmetric=True
    for i in range(numProteinType):
        for j in range(i):
            if (Delta_E[i][j]!=Delta_E[j][i]):
                raiseValueErrorMessage+='''
            "Delta_E" must be a symmetric array!
                '''
                isSymmetric=False
                raiseValueError=True
                break
        if (not(isSymmetric)):
            break
for mass in m:
    if (mass<0):
        raiseValueErrorMessage+='''
        All entries of "m" must be greater than or equal to 0.
        '''
        raiseValueError=True
if (f<=0 or f>1):
    raiseValueErrorMessage+='''
          "f" must be in (0,1]!
    '''
    raiseValueError=True
if (not(type(b)==tuple and len(b)==3)):
    raiseValueErrorMessage+='''
          "b" must be a ternary tuple!
    '''
    raiseValueError=True
else:
    if (int(b[0]-1)!=abs(b[0]-1) or int(b[1]-1)!=abs(b[1]-1) or int(b[2]-1)!=abs(b[2]-1)):
        raiseValueErrorMessage+='''
        All elements of "b" must be positive integers!
        '''
        raiseValueError=True
if (T<=0):
    raiseValueErrorMessage+='''
          "T" must be greater than 0!
    '''
    raiseValueError=True
for n in N:
    if (int(n)!=abs(n)):
        raiseValueErrorMessage+='''
        All entries of "N" must be non-negative integers!
        '''
        raiseValueError=True
    else:
        N=[int(n) for n in N]
if (type(b)==tuple and len(b)==3 and sum(N)>b[0]*b[1]*b[2]):
    raiseValueErrorMessage+='''
          No enough room to place all protein particles!
    '''
    raiseValueError=True
if (N[0]!=0 or N[1]!=0):
    raiseValueErrorMessage+='''
          The number of type "FREE" particles and type "BOUNDARY" particles must be 0!
    '''
    raiseValueError=True
if (R<=0 or k_B<=0):
    raiseValueErrorMessage+='''
          All scientific constants must be positive!
    '''
    raiseValueError=True
if raiseValueError:
    raise ValueError(raiseValueErrorMessage)
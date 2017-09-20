#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
import seaborn as sb

#sb.set(style="white")

#%% 2b Varying dimension for different rhoMax. Parameters
if not os.path.isdir('results'):
    os.mkdir('results')

tolerance = str(1e-9)
numberOfSimulations = str(8)
amplificationFactor = str(2)
maxIterations = str(1e8)
firstH = 0.4

#N = str(25)

rhoMaxVals = ['1', '2.5', '5', '7.5']#, '10'] # For filename
rhoMaxVals2 = [1, 2.5, 5, 7.5]#, 10] # For calculations

#call(["./Allclean"])

#%% 2b Calling cpp 
def RunCpp2b(outfileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit, vectorized):
    counter = 1
    for rhoMax in rhoMaxVals:
        simulationReduction = 0
        N = float(rhoMax)/firstH
        N = int(round(N))
        NVector = [N*(2.)**i for i in xrange(int(numberOfSimulations))]
        for i in NVector:
            if i > NLimit:
                simulationReduction += 1
        numberOfSimulations = int(numberOfSimulations)        
        numberOfSimulations += - simulationReduction
        numberOfSimulations = str(numberOfSimulations)
        N = str(N)
        fileName = outfileName + '%1d' %counter
        print fileName, N
        if NVector[counter - 1] < 800:
            if vectorized:
                call(["./AllrunVectorized", fileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega])
            else:
                call(["./Allrun", fileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega])
            counter += 1
        
outfileName = 'oneElectron'
armadillo = 'false'
electronType = 'oneElectron'
omega = str(1)
rhoMaxVals = ['1', '2.5', '5.0', '7.5']
numberOfSimulations = str(8)
NLimit = 700
vectorized = True
RunCpp2b(outfileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit)

 #%% Plot 1 project 2b   

# Open cpp output
call(["./Allclean"])

oneElectronScalars = {}
for counter in xrange(len(rhoMaxVals2)):
    print counter
    oneElectronScalars[counter+1] = pd.read_table("results/oneElectron%d_scalars.csv" %(counter+1), 
			            delimiter=',')

plt.figure()
legends = []
for key in oneElectronScalars:
    plt.plot(oneElectronScalars[key].h, oneElectronScalars[key].relError)
    legends.append('rhoMax '+rhoMaxVals[key-1])
plt.legend(legends, fontsize = 'large', loc = 0,frameon=False)
plt.title( 'Eigenvalues. Maximum relative errors', fontsize = 'xx-large')
plt.xlabel('h', fontsize = 'xx-large')
plt.ylabel('Max relative error', fontsize = 'xx-large')
plt.grid()
plt.xlim(0.40,0.0)
filename = ('results/oneElectronRelativeErrorEigenvalues.pdf')
plt.savefig(filename)

plt.figure()
legends = []
for key in [3,4]:
    plt.plot(oneElectronScalars[key].h, oneElectronScalars[key].relError)
    legends.append('rhoMax '+rhoMaxVals[key-1])
plt.legend(legends, fontsize = 'large', loc = 0,frameon=False)
plt.title( 'Eigenvalues. Maximum relative errors', fontsize = 'xx-large')
plt.xlabel('h', fontsize = 'xx-large')
plt.ylabel('Max relative error', fontsize = 'xx-large')
plt.grid()
plt.xlim(0.40,0.0)
filename = ('results/oneElectronRelativeErrorEigenvalues2.pdf')
plt.savefig(filename)


#%% Plot 2 project 2b   
plt.figure()
plt.plot(oneElectronScalars[1].N, oneElectronScalars[1].counter)
plt.title( 'Iterations and dimensions', fontsize = 'xx-large')
plt.xlabel('N', fontsize = 'xx-large')
plt.ylabel('Similarity transformations', fontsize = 'xx-large')
plt.grid()

filename = ('results/oneElectronIterationsDimensions.pdf')
plt.savefig(filename)


plt.figure()
plt.plot(oneElectronScalars[1].logN, oneElectronScalars[1].logCounter)
plt.title( 'Iterations and dimensions', fontsize = 'xx-large')
plt.xlabel('log2 N', fontsize = 'xx-large')
plt.ylabel('log2 similarity transformations', fontsize = 'xx-large')
plt.grid()

filename = ('results/oneElectronLogIterationsDimensions.pdf')
plt.savefig(filename)



#%% 2b Comparison Armadillo. Running cpp

outfileName = 'oneElectron'
armadillo = 'true'
electronType = 'oneElectron'
omega = str(1)
rhoMaxVals = ['2.5']
numberOfSimulations = str(8)
NLimit = 700
vectorized = True
RunCpp2b(outfileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit, vectorized)

armadillo = 'false'
vectorized = False
outfileName = 'oneElectronUnvectorized'
RunCpp2b(outfileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit, vectorized)

#%% 2b  Comparison Armadillo plot
# Open cpp output
oneElectronScalarsArmadillo = {}
oneElectronScalarsArmadillo[2] = pd.read_table("results/oneElectron1Armadillo_scalars.csv", 
			            delimiter=',')

oneElectronScalarsUnvectorized = {}
oneElectronScalarsUnvectorized[2] = pd.read_table("results/oneElectronUnvectorized1_scalars.csv", 
			            delimiter=',')

plt.figure()
plt.plot(oneElectronScalarsUnvectorized[2].logN, oneElectronScalarsUnvectorized[2].logTimeUsed, oneElectronScalars[2].logN, oneElectronScalars[2].logTimeUsed, oneElectronScalarsArmadillo[2].logN, oneElectronScalarsArmadillo[2].logTimeUsed)
#plt.title( 'Time used and dimensions\n Jacobi (Vectorized and unvectorized) and armadillo', fontsize = 'xx-large')
plt.legend(['Jacobi unvectorized', 'Jacobi vectorized', 'Armadillo'], fontsize = 'x-large', loc = 0, frameon=False)
plt.xlabel('log N', fontsize = 'xx-large')
plt.ylabel('log time', fontsize = 'xx-large')
plt.grid()
filename = ('results/oneElectronArmadilloLogTimeDimensions.pdf')
plt.savefig(filename)

plt.figure()
plt.plot(oneElectronScalars[2].N, oneElectronScalars[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsArmadillo[2].N, oneElectronScalarsUnvectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed)
plt.legend(['Jacobi vectorized', 'Jacobi unvectorized'], fontsize = 'xx-large', loc = 0)
plt.title( 'Ratio Time Jacobi time Armadillo', fontsize = 'xx-large')
plt.xlabel('N', fontsize = 'xx-large')
plt.ylabel('Time', fontsize = 'xx-large')
plt.grid()

filename = ('results/oneElectronArmadilloTimeDimensions.pdf')
plt.savefig(filename)


#%% project2 d no couloumb interaction
rhoMaxVals = ['5.0']


outfileName = 'twoElectronNoCoulombOmega'
armadillo = 'false'
electronType = 'TwoElectronNoCoulomb'

numberOfSimulations = str(8)
NLimit = 700
vectorized = True

counter = 1
for omega in '0.01', '0.5', '1.0', '5.0':
    fileName = outfileName + '%d' %counter 
    RunCpp2b(fileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit, vectorized)    
    counter += 1

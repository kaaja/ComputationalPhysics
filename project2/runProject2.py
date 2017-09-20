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

rhoMaxVals = ['1', '2.5', '5', '7.5', '10'] # For filename
rhoMaxVals2 = [1, 2.5, 5, 7.5, 10] # For calculations

#call(["./Allclean"])

#%% 2b Calling cpp 
omega = str(1)
counter = 1
for rhoMax in rhoMaxVals:
    N = rhoMaxVals2[counter-1]/firstH
    N = int(round(N))
    N = str(N)
    outfileName = 'oneElectron%1d' %counter
    print outfileName, N
    call(["./AllrunVectorized", outfileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, 'false', 'oneElectron123', omega])
    counter += 1

 #%% Plot 1 project 2b   

# Open cpp output
#call(["./Allclean"])

oneElectronScalars = {}
for counter in xrange(len(rhoMaxVals2)):
    print counter
    oneElectronScalars[counter+1] = pd.read_table("results/oneElectron%d_scalars.csv" %(counter+1), 
			            delimiter=',')

plt.figure()
legends = []
counter = 1
for key in oneElectronScalars:
    plt.plot(oneElectronScalars[key].h, oneElectronScalars[key].relError)
    legends.append('rhoMax '+rhoMaxVals[key-1])
    counter += 1
plt.legend(legends, fontsize = 'large', loc = 'upper right')
plt.title( 'Eigenvalues. Maximum relative errors', fontsize = 'xx-large')
plt.xlabel('h', fontsize = 'xx-large')
plt.ylabel('Max relative error', fontsize = 'xx-large')
plt.grid()
filename = ('results/oneElectronRelativeErrorEigenvalues.pdf')
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
plt.xlabel('log N', fontsize = 'xx-large')
plt.ylabel('log similarity transformations', fontsize = 'xx-large')
plt.grid()

filename = ('results/oneElectronLogIterationsDimensions.pdf')
plt.savefig(filename)



#%% 2b Comparison Armadillo. Running cpp

rhoMax = 2.5
N = rhoMax/firstH
N = int(round(N))
N = str(N)

outfileName = 'oneElectron'
call(["./AllrunVectorized", outfileName, numberOfSimulations,amplificationFactor, N, str(rhoMax), maxIterations, tolerance, 'true', 'oneElecron', omega])

outfileName = 'oneElectronUnvectorized'
call(["./Allrun", outfileName, numberOfSimulations,amplificationFactor, N, str(rhoMax), maxIterations, tolerance, 'false', 'oneElecron', omega])

# Open cpp output
oneElectronScalarsArmadillo = {}
oneElectronScalarsArmadillo[2] = pd.read_table("results/oneElectronArmadillo_scalars.csv", 
			            delimiter=',')

oneElectronScalarsUnvectorized = {}
oneElectronScalarsUnvectorized[2] = pd.read_table("results/oneElectronUnvectorized_scalars.csv", 
			            delimiter=',')

#%% 2b  Comparison Armadillo plot
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
plt.plot(oneElectronScalars[2].N, oneElectronScalars[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalars[2].N, oneElectronScalarsUnvectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed)
plt.legend(['Jacobi vectorized', 'Jacobi unvectorized'], fontsize = 'xx-large', loc = 0)
plt.title( 'Ratio Time Jacobi time Armadillo', fontsize = 'xx-large')
plt.xlabel('N', fontsize = 'xx-large')
plt.ylabel('Time', fontsize = 'xx-large')
plt.grid()

filename = ('results/oneElectronArmadilloTimeDimensions.pdf')
plt.savefig(filename)


#%%  
def plot_logTimes(gaussianTridiagonalScalars, gaussianTridiagonalSymmetricScalars, LUScalars, noLU):
	plt.figure()
	legends = ['Thomas', 'Symmetric']
	plt.plot(gaussianTridiagonalScalars.log_h, gaussianTridiagonalScalars.logTimeUsed)
	plt.hold('on')
	plt.plot(gaussianTridiagonalSymmetricScalars.log_h, gaussianTridiagonalSymmetricScalars.logTimeUsed)
	if not noLU:
	    plt.plot(LUScalars.log_h, LUScalars.logTimeUsed)         
	    legends.append('LU')                 
	plt.legend(legends, fontsize = 'large', loc = 'upper right')
	plt.title('CPU times', fontsize = 'xx-large')
	plt.xlabel('log h', fontsize = 'xx-large')
	plt.ylabel('log time', fontsize = 'xx-large')
	plt.savefig('results/logTimes.pdf')

#%%
def plot_errors(algorithmScalarValues, algorithmNname, amplificationFactor):
	relativeError = algorithmScalarValues.log_rel_error
	plt.figure()
	log_h = algorithmScalarValues.log_h
	plt.plot(log_h, relativeError)
	plt.title('Relative error '+algorithmNname, fontsize = 'xx-large')
	plt.xlabel('log h', fontsize = 'xx-large')
	plt.ylabel('Relative error', fontsize = 'xx-large')
	if int(amplificationFactor) == 2:
	    plt.ylim(-10,-1)        
	filename = ('results/relativeError_'+ algorithmNname + amplificationFactor +'.pdf')
	plt.savefig(filename)


#%% 

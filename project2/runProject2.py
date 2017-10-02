#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict

#sb.set(style="white")

    
#%% 2, run 
def ex2(solverType, tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH, NLimit, vectorized, convergenceLimit, convergenceEigenvalue):
    """
    Calls a function that initializes eigenvalue solvers in cpp.
    Reads .csv files generated from the cpp program and orders the files in dictionaries
    Runs for different omega and rhoMax values
    returns a dict datatype.
    """
    
    electronTypes = ['oneElectron', 'TwoElectronCoulomb', 'TwoElectronNoCoulomb'] #['oneElectron']#
    omegaVals = ['0.01', '0.25', '0.5', '1.0', '5.0']
    rhoMaxVals = ['5','10', '20', '40', '50'] # change back

    
    electronScalars = OrderedDict()
    for electronType in electronTypes:
        electronScalars[electronType]= OrderedDict()
        for omega in omegaVals:
            if (electronType == 'oneElectron' and omega == '1.0') or electronType != 'oneElectron':
                electronScalars[electronType][omega]= OrderedDict()
                counter = 0
                for rhoMax in rhoMaxVals:
                    if solverType == 'armadillo':
                        outfileName = '%sOmega%sRhoMax%sArmadillo' %(electronType,  omega.replace(".", ""), rhoMax.replace(".", ""))
                    elif solverType == 'bisection':
                        outfileName = '%sOmega%sRhoMax%sBisection' %(electronType,  omega.replace(".", ""), rhoMax.replace(".", ""))
                    else:
                        outfileName = '%sOmega%sRhoMax%s' %(electronType,  omega.replace(".", ""), rhoMax.replace(".", ""))
                    print outfileName
                    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
                    electronScalars[electronType][omega][rhoMax] = pd.read_table("results/" + outfileName+ "_scalars.csv", 
        			            delimiter=',')
                    if rhoMax != rhoMaxVals[0]:
                        if np.abs((electronScalars[electronType][omega][rhoMax].lambda1.values[-1] - electronScalars[electronType][omega][rhoMaxVals[counter -1]].lambda1.values[-1])/electronScalars[electronType][omega][rhoMaxVals[counter -1]].lambda1.values[-1]) < convergenceEigenvalue: #' True':
                            break
                    counter += 1
    ex2dPlot(solverType, firstH, electronScalars)
    ex2bplot(firstH, electronScalars['oneElectron'], tolerance, amplificationFactor, maxIterations, numberOfSimulations, convergenceLimit, convergenceEigenvalue)
    
    return electronScalars

#%% 2d, run 
def RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized):
    """
    Compiles and runs cpp program from the command line and
    makes sure the mesh size is not too big.
    """
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
    if vectorized:
        call(["./AllrunVectorized", outfileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit])
    else:
        call(["./Allrun", outfileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit])


#%% 2b, Plotting

def ex2bplot(firstH, electronScalars, tolerance, amplificationFactor, maxIterations, numberOfSimulations, convergenceLimit, convergenceEigenvalue):
    """
    Plots maximum relative error for the one electron case,
    comparison of time used  vs mesh size with different algorithm and/or vectorized and
    log-log plot of time used vs mesh size.
    """
    
    plt.figure()
    legends = []
    for rhoMax in electronScalars['1.0']:
        plt.plot(electronScalars['1.0'][rhoMax].h, electronScalars['1.0'][rhoMax].relError)
        legends.append(r"$\rho_{max}$ = %s" %rhoMax)
    plt.legend(legends, loc = 0)
    plt.title( 'Eigenvalues. Maximum relative errors')
    plt.xlabel('h')
    plt.ylabel('Max relative error')
    plt.xlim(0.40,0.0)
    filename = ('results/oneElectronRelativeErrorEigenvalues.pdf')
    plt.savefig(filename)
    plt.close()

    # Open cpp output
    firstH = 0.4
    electronType = 'oneElectron'
    omega = '1.0'
    rhoMax = '5'
    NLimit = 500
    
    outfileName = 'oneElectronArmadilloTimeComparison'
    solverType = 'armadillo'
    vectorized = True
    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronUnvectorizedTimeComparison'
    solverType = 'jacobi'
    vectorized = False
    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronVectorizedTimeComparison'
    solverType = 'jacobi'
    vectorized = True
    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronBisectionVectorizedTimeComparison'
    solverType = 'bisection'
    vectorized = True
    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronLanczosVectorizedTimeComparison'
    solverType = 'lanczosArmadillo'
    vectorized = True
    RunCpp(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, solverType, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    
    oneElectronScalarsUnvectorized = {}
    oneElectronScalarsUnvectorized[2] = pd.read_table("results/oneElectronUnvectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronScalarsArmadillo = {}
    oneElectronScalarsArmadillo[2] = pd.read_table("results/oneElectronArmadilloTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronScalarsVectorized = {}
    oneElectronScalarsVectorized[2] = pd.read_table("results/oneElectronVectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronBisectionScalarsVectorized = {}
    oneElectronBisectionScalarsVectorized[2] = pd.read_table("results/oneElectronBisectionVectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronLanczosScalarsVectorized = {}
    oneElectronLanczosScalarsVectorized[2] = pd.read_table("results/oneElectronLanczosVectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    
    plt.figure()
    plt.plot(oneElectronLanczosScalarsVectorized[2].logN, oneElectronLanczosScalarsVectorized[2].logTimeUsed, oneElectronBisectionScalarsVectorized[2].logN, oneElectronBisectionScalarsVectorized[2].logTimeUsed, oneElectronScalarsUnvectorized[2].logN, oneElectronScalarsUnvectorized[2].logTimeUsed, oneElectronScalarsVectorized[2].logN, oneElectronScalarsVectorized[2].logTimeUsed, oneElectronScalarsArmadillo[2].logN, oneElectronScalarsArmadillo[2].logTimeUsed)
    #plt.title( 'Time used and dimensions\n Jacobi (Vectorized and unvectorized) and armadillo', fontsize = 'xx-large')
    plt.legend(['Lanczos vectorized','Bisection vectorized', 'Jacobi unvectorized', 'Jacobi vectorized', 'Armadillo'], loc = 0)
    plt.xlabel('log N')
    plt.ylabel('log time')
    filename = ('results/oneElectronComparisonLogTimeDimensions.pdf')
    plt.savefig(filename)
    
    plt.figure()
    plt.plot(oneElectronLanczosScalarsVectorized[2].N, oneElectronLanczosScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronBisectionScalarsVectorized[2].N, oneElectronBisectionScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsVectorized[2].N, oneElectronScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsArmadillo[2].N, oneElectronScalarsUnvectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed)
    plt.legend(['Lanczos vectorized', 'Bisection vectorized', 'Jacobi vectorized', 'Jacobi unvectorized'], loc = 0)
    plt.title( 'Time relative to armadillo time')
    plt.xlabel('N')
    plt.ylabel('Time')
    filename = ('results/oneElectronArmadilloTimeDimensions.pdf')
    plt.savefig(filename)
    plt.close()
    

#%% 2d, Plots 
def ex2dPlot(solverType, firstH, electronScalars):        
    """
    Function makes 2 kinds of plots:
    1) Plots minimum eigenvalues vs mesh size for a specific frequency, and comparing with increasing 
    rhoMax values.
    
    2) Plots minimum eigenvalue vs mesh size and comparing for different frequency for a suitable 
    rhoMax value.
    
    """
    
    
    for electronType in electronScalars:
        fig2, ax2 = plt.subplots()
        legends2 = []
        for omega in electronScalars[electronType]:
            fig, ax = plt.subplots()
            legends = []
            for rhoMax in electronScalars[electronType][omega]:
                if omega == '0.25' and electronType == 'TwoElectronCoulomb':
                    labels = 'Relative error'
                    paperValue = 1.25
                    relativeError = np.abs((electronScalars[electronType][omega][rhoMax].lambda1 - paperValue)/paperValue)
                    ax.plot(electronScalars[electronType][omega][rhoMax].h, relativeError)
                else:
                    labels = 'Minimum ' + r'$\lambda$'
                    ax.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                legends.append(r"$\rho_{max}$ = %s" %rhoMax)
            ax.legend(legends, fontsize = 'large', loc = 0,frameon=False)
            ax.set_title( "Minimum $\lambda$, %s $\omega$ = %s" %(electronType, omega))
            ax.set_xlabel('h')
            ax.set_ylabel(labels)
            ax.grid(True)
            ax.set_xlim(float(firstH),0.0)
            #fig.show()
            outfileName = '%sOmega%sRhoMaxComparison' %(electronType,  omega.replace(".", ""))
            filename = ('results/' + outfileName + '.pdf')
            fig.savefig(filename) 
            plt.close()
            
            # Plots different omega's same figure
            ax2.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
            ax2.hold('on')
            legends2.append(r"$\omega$ = %s" %omega)

            
        ax2.legend(legends2,loc = 0)
        ax2.set_title( "Minimum $\lambda$, %s" %electronType)
        ax2.set_xlabel('h')
        ax2.set_ylabel(labels)
        ax2.set_xlim(float(firstH),0.0)
        #fig.show()
        outfileName = '%sOmegaComparison' %(electronType)
        filename = ('results/' + outfileName + '.pdf')
        fig2.savefig(filename)
        plt.close()

#%% Command line options 

if __name__ == "__main__":
	try:
		import seaborn
		seaborn.set(style="white", context="notebook", font_scale=1.5,
		            rc={"axes.grid": True, "legend.frameon": False,
		                "lines.markeredgewidth": 1.4, "lines.markersize": 10})
	except ImportError:
		print("Could not import seaborn for plot styling. Try")
		print("\n    conda install seaborn\n\nor")
		print("\n    pip install seaborn\n")
		
	tolerance = str(1e-6) # Kja: testing with low. CHange back if problems
	numberOfSimulations = str(16)
	amplificationFactor = str(2.0)
	maxIterations = str(1e8)
	firstH = 0.4
	NLimit = 2000
	vectorized = True
	convergenceLimit = '.0011' # for N
	convergenceEigenvalue = .001 # For use in python, deciding wheter new rhoMax Needed
	
	
	
	parser = argparse.ArgumentParser(description="starts a c++ program solving eigenvalue problems, reads and  plots.")
	parser.add_argument("task", type=str, default='2b', help="choose task to solve. 2b, 2b Armadillo or 2d")
	parser.add_argument("algorithm", type=str, default='armadillo', help="choose algorithm to solve eigenvalue problem. jacobi, armadillo, bisection or lanczosArmadillo")
	args = parser.parse_args()

	if not os.path.isdir('results'):
		os.mkdir('results')
		
	if args.task == '2':
		electronScalars = ex2(args.algorithm, tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH, NLimit, vectorized, convergenceLimit, convergenceEigenvalue)


    

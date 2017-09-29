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
    
    electronTypes = ['oneElectron', 'TwoElectronCoulomb', 'TwoElectronNoCoulomb']
    omegaVals = ['0.01', '0.25', '0.5', '1.0', '5.0']#'0.25']#, '0.5', '1.0', '5.0']
    rhoMaxVals = ['5','10', '20', '40', '50']

    
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
    
    
    
    plt.figure()
    plt.plot(oneElectronBisectionScalarsVectorized[2].logN, oneElectronBisectionScalarsVectorized[2].logTimeUsed, oneElectronScalarsUnvectorized[2].logN, oneElectronScalarsUnvectorized[2].logTimeUsed, oneElectronScalarsVectorized[2].logN, oneElectronScalarsVectorized[2].logTimeUsed, oneElectronScalarsArmadillo[2].logN, oneElectronScalarsArmadillo[2].logTimeUsed)
    #plt.title( 'Time used and dimensions\n Jacobi (Vectorized and unvectorized) and armadillo', fontsize = 'xx-large')
    plt.legend(['Bisection vectorized', 'Jacobi unvectorized', 'Jacobi vectorized', 'Armadillo'], loc = 0)
    plt.xlabel('log N')
    plt.ylabel('log time')
    filename = ('results/oneElectronComparisonLogTimeDimensions.pdf')
    plt.savefig(filename)
    
    plt.figure()
    plt.plot(oneElectronBisectionScalarsVectorized[2].N, oneElectronBisectionScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsVectorized[2].N, oneElectronScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsArmadillo[2].N, oneElectronScalarsUnvectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed)
    plt.legend(['Bisection vectorized', 'Jacobi vectorized', 'Jacobi unvectorized'], loc = 0)
    plt.title( 'Time relative to armadillo time')
    plt.xlabel('N')
    plt.ylabel('Time')
    filename = ('results/oneElectronArmadilloTimeDimensions.pdf')
    plt.savefig(filename)
    plt.close()
    
#ex2bplot('true', '0.4', electronScalars['oneElectron'])

#%% 2d, Plots 
def ex2dPlot(solverType, firstH, electronScalars):        
    """
    Function makes 3 kinds of plots:
    1) Plots minimum eigenvalues vs mesh size for a specific frequency, and comparing with increasing 
    rhoMax values.
    
    2) Plots minimum eigenvalue vs mesh size and comparing for different frequency for a suitable 
    rhoMax value.
    
    3) Plots minimum eigenvalue vs mesh size and comparing one elctron, two electron with no interaction and
    two electron with interaction for given frequency and suitable rhoMax value.
    """
    
    fig3, ax3 = plt.subplots() # For different electron type in same plot
    ax3.hold('on')
    
    fig4, ax4 = plt.subplots() # For different electron type in same plot
    ax4.hold('on')
    
    fig5, ax5 = plt.subplots() # For different electron type in same plot
    ax5.hold('on')
    
    fig6, ax6 = plt.subplots() # For different electron type in same plot
    ax6.hold('on')
    
    fig7, ax7 = plt.subplots() # For different electron type in same plot
    ax7.hold('on')
    
    legends3 = []
    
    for electronType in electronScalars:
        fig2, ax2 = plt.subplots()
        legends2 = []
        legends3.append(electronType)
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
            
            if electronType != 'oneElectron':
                if omega == '0.01': 
                    ax3.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                    ax3.set_title(r"$\omega$ = %s" %omega)
                if omega == '0.25': 
                    ax4.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                    ax4.set_title(r"$\omega$ = %s" %omega)
                if omega == '0.5': 
                    ax5.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                    ax5.set_title(r"$\omega$ = %s" %omega)
                if omega == '5.0': 
                    ax7.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                    ax7.set_title(r"$\omega$ = %s" %omega)
            else:
                if omega == '1.0': 
                    ax6.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                    ax6.set_title(r"$\omega$ = %s" %omega)

            
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
        
    ax3.legend(legends3[1:3], loc = 0)
    ax3.set_xlabel('h')
    ax3.set_ylabel(labels)
    ax3.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = 'electronComparisonOmega001' 
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    plt.close()

    ax4.legend(legends3[1:3], loc = 0)
    ax4.set_xlabel('h')
    ax4.set_ylabel(labels)
    ax4.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = 'electronComparisonOmega025' 
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    plt.close()

    ax5.legend(legends3[1:3],  loc = 0)
    ax5.set_xlabel('h')
    ax5.set_ylabel(labels)
    ax5.grid(True)
    ax5.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = 'electronComparisonOmega05'
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    plt.close()
    
    ax6.legend(legends3, loc = 0)
    ax6.set_xlabel('h')
    ax6.set_ylabel(labels)
    ax6.grid(True)
    ax6.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = 'electronComparisonOmega50'
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    plt.close()
    
    ax7.legend(legends3[1:3], loc = 0)
    ax7.set_xlabel('h')
    ax7.set_ylabel(labels)
    ax7.grid(True)
    ax7.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = 'electronComparisonOmega10'
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
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
		
	solverType = 'armadillo'
	tolerance = str(1e-6) # Kja: testing with low. CHange back if problems
	numberOfSimulations = str(16)
	amplificationFactor = str(2.0)
	maxIterations = str(1e8)
	firstH = 0.4
	NLimit = 2000
	vectorized = True
	convergenceLimit = '.001'
	convergenceEigenvalue = .001 # For use in python, deciding wheter new rhoMax Needed
	
	
	
	parser = argparse.ArgumentParser(description="starts a c++ program solving eigenvalue problems, reads and  plots.")
	parser.add_argument("task", type=str, default='2b', help="choose task to solve. 2b, 2b Armadillo or 2d")
	args = parser.parse_args()

	if not os.path.isdir('results'):
		os.mkdir('results')

	if  args.task == '2b':
	   ex2b(tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH)
	elif args.task == '2bArmadillo':
		ex2barmadillo(tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH)
	elif args.task == '2':
		electronScalars = ex2(solverType, tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH, NLimit, vectorized, convergenceLimit, convergenceEigenvalue)


    

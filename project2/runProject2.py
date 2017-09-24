#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
import seaborn as sb
from collections import OrderedDict

#sb.set(style="white")


#%% Function calling cpp
def RunCpp2b(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega,convergenceLimit, NLimit, vectorized):
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
                call(["./AllrunVectorized", fileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit])
            else:
                call(["./Allrun", fileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit])
            counter += 1


#%%
def ex2b():
    outfileName = 'oneElectron'
    armadillo = 'false'
    electronType = 'oneElectron'
    omega = str(1)
    rhoMaxVals = ['1', '2.5', '5.0', '7.5']
    numberOfSimulations = str(8)
    NLimit = 700
    vectorized = True
    RunCpp2b(outfileName, numberOfSimulations,amplificationFactor, rhoMaxVals, maxIterations, tolerance, armadillo, electronType, omega, NLimit)

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
    
#%% 2b, Armadillo comparison 

def ex2barmadillo():
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
    
#%% 2b, Plotting

def ex2bplot(armadillo, firstH, electronScalars):
    
    plt.figure()
    legends = []
    for rhoMax in electronScalars['1.0']:
        plt.plot(electronScalars['1.0'][rhoMax].h, electronScalars['1.0'][rhoMax].relError)
        legends.append('rhoMax %s' %rhoMax)
    plt.legend(legends, fontsize = 'large', loc = 0,frameon=False)
    plt.title( 'Eigenvalues. Maximum relative errors', fontsize = 'xx-large')
    plt.xlabel('h', fontsize = 'xx-large')
    plt.ylabel('Max relative error', fontsize = 'xx-large')
    plt.grid()
    plt.xlim(0.40,0.0)
    filename = ('results/oneElectronRelativeErrorEigenvalues.pdf')
    plt.savefig(filename)
    
    
#    plt.figure()
#    rhoMax = '10'
#    plt.plot(electronScalars['1.0'][rhoMax].N, electronScalars['1.0'][rhoMax].counter)
#    plt.title( 'Iterations and dimensions', fontsize = 'xx-large')
#    plt.xlabel('N', fontsize = 'xx-large')
#    plt.ylabel('Similarity transformations', fontsize = 'xx-large')
#    plt.grid()
#    filename = ('results/oneElectronIterationsDimensions.pdf')
#    plt.savefig(filename)
#     
#    plt.figure()
#    plt.plot(electronScalars['1.0'][rhoMax].logN, electronScalars['1.0'][rhoMax].logCounter)
#    plt.title( 'Iterations and dimensions', fontsize = 'xx-large')
#    plt.xlabel('log2 N', fontsize = 'xx-large')
#    plt.ylabel('log2 similarity transformations', fontsize = 'xx-large')
#    plt.grid()
#    filename = ('results/oneElectronLogIterationsDimensions.pdf')
#    plt.savefig(filename)

    # Open cpp output
    tolerance = str(1e-9)
    amplificationFactor = str(2)
    maxIterations = str(1e8)
    firstH = 0.4
    electronType = 'oneElectron'
    omega = '1.0'
    rhoMax = '5'
    numberOfSimulations = str(15)
    NLimit = 2000
    convergenceLimit = '.0001'
    convergenceEigenvalue = .000001
    
    outfileName = 'oneElectronArmadilloTimeComparison'
    armadillo = 'true'
    vectorized = True
    RunCpp2d(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronUnvectorizedTimeComparison'
    armadillo = 'false'
    vectorized = False
    RunCpp2d(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    outfileName = 'oneElectronVectorizedTimeComparison'
    armadillo = 'false'
    vectorized = True
    RunCpp2d(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit, NLimit, vectorized)
    
    oneElectronScalarsUnvectorized = {}
    oneElectronScalarsUnvectorized[2] = pd.read_table("results/oneElectronUnvectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronScalarsArmadillo = {}
    oneElectronScalarsArmadillo[2] = pd.read_table("results/oneElectronArmadilloTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    oneElectronScalarsVectorized = {}
    oneElectronScalarsVectorized[2] = pd.read_table("results/oneElectronVectorizedTimeComparison_scalars.csv", 
    			            delimiter=',')
    
    
    plt.figure()
    plt.plot(oneElectronScalarsUnvectorized[2].logN, oneElectronScalarsUnvectorized[2].logTimeUsed, oneElectronScalarsVectorized[2].logN, oneElectronScalarsVectorized[2].logTimeUsed, oneElectronScalarsArmadillo[2].logN, oneElectronScalarsArmadillo[2].logTimeUsed)
    #plt.title( 'Time used and dimensions\n Jacobi (Vectorized and unvectorized) and armadillo', fontsize = 'xx-large')
    plt.legend(['Jacobi unvectorized', 'Jacobi vectorized', 'Armadillo'], fontsize = 'x-large', loc = 0, frameon=False)
    plt.xlabel('log N', fontsize = 'xx-large')
    plt.ylabel('log time', fontsize = 'xx-large')
    plt.grid()
    filename = ('results/oneElectronComparisonLogTimeDimensions.pdf')
    plt.savefig(filename)
    
    plt.figure()
    plt.plot(oneElectronScalarsVectorized[2].N, oneElectronScalarsVectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed, oneElectronScalarsArmadillo[2].N, oneElectronScalarsUnvectorized[2].timeUsed/oneElectronScalarsArmadillo[2].timeUsed)
    plt.legend(['Jacobi vectorized', 'Jacobi unvectorized'], fontsize = 'xx-large', loc = 0)
    plt.title( 'Ratio Time Jacobi time Armadillo', fontsize = 'xx-large')
    plt.xlabel('N', fontsize = 'xx-large')
    plt.ylabel('Time', fontsize = 'xx-large')
    plt.grid()
    filename = ('results/oneElectronArmadilloTimeDimensions.pdf')
    plt.savefig(filename)
    
ex2bplot('true', '0.4', electronScalars['oneElectron'])
    
#%% 2d, run 
def ex2d():
    armadillo = 'true'
    tolerance = str(1e-9)
    numberOfSimulations = str(8)
    amplificationFactor = str(2)
    maxIterations = str(1e8)
    firstH = 0.4
    electronTypes = ['oneElectron', 'TwoElectronCoulomb', 'TwoElectronNoCoulomb']
    omegaVals = ['0.01', '1.0']#'0.25']#, '0.5', '1.0', '5.0']
    rhoMaxVals = ['5','10']#, '15', '20', '40','50'] 
    numberOfSimulations = str(15)
    NLimit = 2000
    vectorized = True
    convergenceLimit = '.0001'
    convergenceEigenvalue = .000001
    
    electronScalars = OrderedDict()
    for electronType in electronTypes:
        electronScalars[electronType]= OrderedDict()
        for omega in omegaVals:
            if (electronType == 'oneElectron' and omega == '1.0') or electronType != 'oneElectron':
                electronScalars[electronType][omega]= OrderedDict()
                counter = 0
                for rhoMax in rhoMaxVals:
                    if armadillo == 'true':
                        outfileName = '%sOmega%sRhoMax%sArmadillo' %(electronType,  omega.replace(".", ""), rhoMax.replace(".", ""))
                    else:
                        outfileName = '%sOmega%sRhoMax%s' %(electronType,  omega.replace(".", ""), rhoMax.replace(".", ""))
                    print outfileName
                    RunCpp2d(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit, NLimit, vectorized)
                    electronScalars[electronType][omega][rhoMax] = pd.read_table("results/" + outfileName+ "_scalars.csv", 
        			            delimiter=',')
                    if rhoMax != rhoMaxVals[0]:
                        if np.abs((electronScalars[electronType][omega][rhoMax].lambda1.values[-1] - electronScalars[electronType][omega][rhoMaxVals[counter -1]].lambda1.values[-1])/electronScalars[electronType][omega][rhoMaxVals[counter -1]].lambda1.values[-1]) < convergenceEigenvalue: #' True':
                            break
                    counter += 1
    ex2dPlot(armadillo, firstH, electronScalars)
    
    return electronScalars

#%% 2d, run 
def RunCpp2d(outfileName, firstH, numberOfSimulations,amplificationFactor, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit, NLimit, vectorized):
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
        call(["./AllrunVectorized", outfileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit])
    else:
        call(["./Allrun", outfileName, numberOfSimulations,amplificationFactor, N, rhoMax, maxIterations, tolerance, armadillo, electronType, omega, convergenceLimit])

#%% 2d, Plots 
def ex2dPlot(armadillo, firstH, electronScalars):        
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
                    labels = 'Minimum eigenvalue'
                    ax.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                legends.append('rhoMax %s'%rhoMax)
            ax.legend(legends, fontsize = 'large', loc = 0,frameon=False)
            ax.set_title( 'Minimum eigenvalues. \n %s Omega %s' %(electronType, omega), fontsize = 'xx-large')
            ax.set_xlabel('h', fontsize = 'xx-large')
            ax.set_ylabel(labels, fontsize = 'xx-large')
            ax.grid(True)
            ax.set_xlim(float(firstH),0.0)
            #fig.show()
            outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
            filename = ('results/' + outfileName + '.pdf')
            fig.savefig(filename) 
            
            # Plots different omega's same figure
            ax2.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
            ax2.hold('on')
            legends2.append('Omega = %s' %omega)
            
            if omega == '0.01': 
                ax3.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                ax3.set_title('Omega = %s' %omega, fontsize = 'xx-large')
            if omega == '0.25': 
                ax4.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                ax4.set_title('Omega = %s' %omega, fontsize = 'xx-large')
            if omega == '0.5': 
                ax5.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                ax5.set_title('Omega = %s' %omega, fontsize = 'xx-large')
            if omega == '1.0': 
                ax6.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                ax6.set_title('Omega = %s' %omega, fontsize = 'xx-large')
            if omega == '5.0': 
                ax7.plot(electronScalars[electronType][omega][rhoMax].h, electronScalars[electronType][omega][rhoMax].lambda1)
                ax7.set_title('Omega = %s' %omega, fontsize = 'xx-large')
            
            
        ax2.legend(legends2, fontsize = 'large', loc = 0,frameon=False)
        ax2.set_title( 'Minimum eigenvalue. \n %s ' %electronType, fontsize = 'xx-large')
        ax2.set_xlabel('h', fontsize = 'xx-large')
        ax2.set_ylabel(labels, fontsize = 'xx-large')
        ax2.grid(True)
        ax2.set_xlim(float(firstH),0.0)
        #fig.show()
        outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
        filename = ('results/' + outfileName + '.pdf')
        fig2.savefig(filename)
    
        
    ax3.legend(legends3, fontsize = 'large', loc = 0,frameon=False)
    ax3.set_xlabel('h', fontsize = 'xx-large')
    ax3.set_ylabel(labels, fontsize = 'xx-large')
    ax3.grid(True)
    ax3.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)

    ax4.legend(legends3, fontsize = 'large', loc = 0,frameon=False)
    ax4.set_xlabel('h', fontsize = 'xx-large')
    ax4.set_ylabel(labels, fontsize = 'xx-large')
    ax4.grid(True)
    ax4.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)

    ax5.legend(legends3, fontsize = 'large', loc = 0,frameon=False)
    ax5.set_xlabel('h', fontsize = 'xx-large')
    ax5.set_ylabel(labels, fontsize = 'xx-large')
    ax5.grid(True)
    ax5.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    
    ax6.legend(legends3, fontsize = 'large', loc = 0,frameon=False)
    ax6.set_xlabel('h', fontsize = 'xx-large')
    ax6.set_ylabel(labels, fontsize = 'xx-large')
    ax6.grid(True)
    ax6.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)
    
    ax7.legend(legends3, fontsize = 'large', loc = 0,frameon=False)
    ax7.set_xlabel('h', fontsize = 'xx-large')
    ax7.set_ylabel(labels, fontsize = 'xx-large')
    ax7.grid(True)
    ax7.set_xlim(float(firstH),0.0)
    #fig.show()
    outfileName = '%sOmega%s' %(electronType,  omega.replace(".", ""))
    filename = ('results/' + outfileName + '.pdf')
    fig3.savefig(filename)


    #ex2dPlot('true', '0.4', electronScalars)

#%% Command line options 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="starts a c++ program solving eigenvalue problems, reads and  plots.")
    parser.add_argument("task", type=str, default='2b', help="choose task to solve. 2b, sbArmadillo or 2d")
    args = parser.parse_args()
    
    if not os.path.isdir('results'):
        os.mkdir('results')
    
    if  args.task == '2b':
       ex2b(tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH)
    elif args.task == '2bArmadillo':
        ex2barmadillo(tolerance, numberOfSimulations, amplificationFactor, maxIterations, firstH)
    elif args.task == '2d':
        electronScalars = ex2d()
    
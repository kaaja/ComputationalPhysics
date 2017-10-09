#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict

try:
    import seaborn
    seaborn.set(style="white", context="notebook", font_scale=1.5,
		            rc={"axes.grid": True, "legend.frameon": False,
		                "lines.markeredgewidth": 1.4, "lines.markersize": 10})
except ImportError:
	print("Could not import seaborn for plot styling. Try")
	print("\n    conda install seaborn\n\nor")
	print("\n    pip install seaborn\n")

if not os.path.isdir('results'):
	os.mkdir('results')

#%% Run 

def runCpp(outfileName, finalTime, N, solverType):
    """
    Compiles and runs cpp program from the command line and
    makes sure the mesh size is not too big.
    """
    N = str(N)
    finalTime = str(finalTime)
    call(["./AllrunVectorized", outfileName, finalTime, N, solverType])
    
#%% 2, run 

def sunEarth():
    """

    """
    sunEarth      = OrderedDict()
    supNormValues = []#OrderedDict()
    supNormAngularMomentum = []
    
    outfileName = 'sunEarth'
    solverType = 'forwardEuler'
    
    finalTimes = [10**i for i in xrange(0,4)]
    #Ns = [10**i for i in xrange(3,8)]   
    dts = [10.**(-i) for i in xrange(1, 5)]
    
    
    # Plots. Positions vs time
    fig, ax = plt.subplots(2, sharex=True)
    fig.hold('on')
    ax[0].set_title(outfileName + ' ' + solverType)
    ax[1].set_xlabel('t [Au]')
    ax[0].set_ylabel('x [Au]')
    ax[1].set_ylabel('y [Au]')
    plt.ylim(-1.1, 1.1)    
    legends = []
    
    # Energy plots
    fig2, ax2 = plt.subplots()
    fig2.hold('on')
    ax2.set_title(outfileName + ' ' + solverType + 'Energy')
    ax2.set_xlabel('t [Au]')
    ax2.set_ylabel('E ')
    ax2.set_ylim(0.7, 1.05)    
    

    
    # Angular momentum plots
    fig3, ax3 = plt.subplots()
    fig3.hold('on')
    ax3.set_title(outfileName + ' ' + solverType + 'Angular momentum')
    ax3.set_xlabel('t [Au]')
    ax3.set_ylabel('Angular momentum ')
    ax3.set_ylim(0.99, 1.05)    
    
    # supNorm ENergy plots
    fig4, ax4 = plt.subplots()
    fig4.hold('on')
    ax4.set_title(outfileName + ' ' + solverType + 'sup Norm')
    ax4.set_xlabel("$\Delta t$")
    ax4.set_ylabel('sup norm Energy')
    #ax2.set_ylim(0.)

    # supNorm Angular Momentum
    fig5, ax5 = plt.subplots()
    fig5.hold('on')
    ax5.set_title(outfileName + ' ' + solverType + ' sup-norm Angular momentum')
    ax5.set_xlabel("$\Delta t$")
    ax5.set_ylabel('sup-norm angular momentum')
    #ax2.set_ylim(0.)
    
    for finalTime in finalTimes:
        sunEarth['FinalTime %f' %finalTime] = {}
        for dt in dts:
            N = finalTime/dt
            print 'N = %d, final time = %.2g' %(N, finalTime)  
            outfileName2 = outfileName + 'finalTime%s' %str(finalTime).replace(".", "") + 'N%s' %str(N).replace(".", "")
            runCpp(outfileName2, finalTime, N, solverType)
            sunEarth['FinalTime %f' %finalTime]['N %f' %N] = pd.read_table("results/" + outfileName2 + ".csv", 
            			            delimiter=',')
            plotSunEarth(sunEarth['FinalTime %f' %finalTime]['N %f' %N], outfileName2, solverType, N, finalTime)
            if finalTime == finalTimes[1]:
                #ax2.loglog(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy)
                ax2.plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, (sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy)/(sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy[0] + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy[0]))
                ax3.plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, (sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum )/(sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum[0]))
                ax[0].plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, sunEarth['FinalTime %f' %finalTime]['N %f' %N].x)
                ax[1].plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, sunEarth['FinalTime %f' %finalTime]['N %f' %N].y)
                #fig.savefig('results/' + outfileName + 'times.png') 
                #plt.show()
                #plt.close()
                legends.append('dt %.2g' %dt)
            
            if finalTime == finalTimes[-1]:
                totalEnergy = (sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy)
                totalEnergy = totalEnergy.values 
                angularMomentum = (sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum).values
                supNormValues.append(supNorm(totalEnergy))
                supNormAngularMomentum.append(supNorm(angularMomentum))
    
    ax[0].legend(legends)
    ax[1].legend(legends)
    fig.savefig('results/' + outfileName + 'Times.png') 
    
    ax2.legend(legends)
    fig2.savefig('results/'+ outfileName + 'Energy.png')
    
    ax3.legend(legends)
    fig3.savefig('results/'+ outfileName + 'AngularMomentum.png')
    
    ax4.loglog(dts, supNormValues)
    fig4.savefig('results/'+ outfileName + 'supNorm.png')
    
    ax5.loglog(dts, supNormAngularMomentum)
    fig5.savefig('results/'+ outfileName + 'supNormAngularMomentum.png')
    
    return sunEarth, supNormValues, supNormAngularMomentum

def supNorm(totalEnergy):
    totalEnergy = totalEnergy[1:]/totalEnergy[0]-1
    supNorm = np.max(abs(totalEnergy))
    return supNorm


def plotSunEarth(sunEarth, outfileName, solverType, N, finalTime):
    dt = finalTime/N
    fig, ax = plt.subplots()
    ax.plot(sunEarth.x, sunEarth.y)
    ax.set_title(outfileName + ' ' + solverType + '\n N %d, dt %.2g' %(N, dt))
    ax.set_xlabel('x [Au]')
    ax.set_ylabel('y [Au]')
    ax.set_xlim(-2., 2.)
    ax.set_ylim(-2., 2.)
    fig.savefig('results/' + outfileName + '.png') 
    #plt.show()
    plt.close()
    

def plotSunEarthTimes(sunEarth, outfileName, solverType, N, finalTime):
    dt = finalTime/N
    fig, ax = plt.subplots(2, sharex=True)
    fig.hold('on')
    ax[0].plot(sunEarth.time, sunEarth.x)
    ax[1].plot(sunEarth.time, sunEarth.y)
    fig.title(outfileName + ' ' + solverType + '\n N %d, dt %.2g' %(N, dt))
    ax[0].set_xlabel('t [Au]')
    ax[0].set_ylabel('x [Au]')
    ax[1].set_ylabel('y [Au]')
    #ax.set_xlim(-2., 2.)
    fig.ylim(-1.1, 1.1)
    #fig.savefig('results/' + outfileName + 'times.png') 
    #plt.show()
    #plt.close()


#%% 2, run     
sunearth, supNormValues, supNormAngularMomentum = sunEarth()

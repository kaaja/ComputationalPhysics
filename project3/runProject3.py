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

def runCpp(outfileName, finalTime, N, solverType, initialVy):
    """
    Compiles and runs cpp program from the command line and
    makes sure the mesh size is not too big.
    """
    N = str(N)
    finalTime = str(finalTime)
    initialVy = str(initialVy)
    call(["./AllrunVectorized", outfileName, finalTime, N, solverType, initialVy])
    
#%% 2, run 

def sunEarth():
    """

    """
    initialVy = 2*np.pi
    sunEarth      = OrderedDict()
    supNormValues = []#OrderedDict()
    supNormAngularMomentum = []
    timesUsed = []
    
    outfileName = 'sunEarth'
    solverType = 'VelocityVerlet'
    
    finalTimes = [10**i for i in xrange(0,4)]
    #Ns = [10**i for i in xrange(3,8)]   
    dts = [10.**(-i) for i in xrange(1, 5)]
    
    Ns = np.asarray(finalTimes)/np.asarray(dts)
    
    
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
    
    # Times
    fig6, ax6 = plt.subplots()
    fig6.hold('on')
    ax6.set_title(outfileName + ' ' + solverType + ' log times')
    ax6.set_xlabel("log N")
    ax6.set_ylabel('log time')
    
    for finalTime in finalTimes:
        sunEarth['FinalTime %f' %finalTime] = {}
        for dt in dts:
            N = finalTime/dt
            print 'N = %d, final time = %.2g' %(N, finalTime)  
            outfileName2 = outfileName + 'finalTime%s' %str(finalTime).replace(".", "") + 'N%s' %str(N).replace(".", "")
            runCpp(outfileName2, finalTime, N, solverType, initialVy)
            sunEarth['FinalTime %f' %finalTime]['N %f' %N] = pd.read_table("results/" + outfileName2 + ".csv", 
            			            delimiter=',')
            plotSunEarth(sunEarth['FinalTime %f' %finalTime]['N %f' %N], outfileName2, solverType, N, finalTime)
            if finalTime == finalTimes[1]:
                #ax2.loglog(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time, sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy)
                ax2.plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time[:-1], (sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy[:-1] + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy[:-1])/(sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy[0] + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy[0]))
                ax3.plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time[:-1], (sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum[:-1] )/(sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum[0]))
                ax[0].plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time[:-1], sunEarth['FinalTime %f' %finalTime]['N %f' %N].x[:-1])
                ax[1].plot(sunEarth['FinalTime %f' %finalTime]['N %f' %N].time[:-1], sunEarth['FinalTime %f' %finalTime]['N %f' %N].y[:-1])
                #fig.savefig('results/' + outfileName + 'times.png') 
                #plt.show()
                #plt.close()
                legends.append('dt %.2g' %dt)
            
            if finalTime == finalTimes[-1]:
                totalEnergy = (sunEarth['FinalTime %f' %finalTime]['N %f' %N].kineticEnergy[:-1] + sunEarth['FinalTime %f' %finalTime]['N %f' %N].potentialEnergy[:-1])
                totalEnergy = totalEnergy.values 
                angularMomentum = (sunEarth['FinalTime %f' %finalTime]['N %f' %N].angularMomentum[:-1]).values
                supNormValues.append(supNorm(totalEnergy))
                supNormAngularMomentum.append(supNorm(angularMomentum))
                timesUsed.append(sunEarth['FinalTime %f' %finalTime]['N %f' %N].timeUsed.values[-1])
    
    ax[0].legend(legends)
    ax[1].legend(legends)
    fig.savefig('results/' + outfileName + 'Times' + solverType + '.png') 
    
    ax2.legend(legends)
    fig2.savefig('results/'+ outfileName + 'Energy' + solverType + '.png')
    
    ax3.legend(legends)
    fig3.savefig('results/'+ outfileName + 'AngularMomentum' + solverType + '.png')
    
    ax4.loglog(dts, supNormValues)
    fig4.savefig('results/'+ outfileName + 'supNorm' + solverType + '.png')
    
    ax5.loglog(dts, supNormAngularMomentum)
    fig5.savefig('results/'+ outfileName + 'supNormAngularMomentum' + solverType + '.png')
    
    ax6.loglog(Ns, timesUsed)
    fig6.savefig('results/'+ outfileName + 'logTimeUsed' + solverType + '.png')
    
    plt.close()
    
    return sunEarth, supNormValues, supNormAngularMomentum

def supNorm(totalEnergy):
    totalEnergy = totalEnergy[1:]/totalEnergy[0]-1
    supNorm = np.max(abs(totalEnergy))
    return supNorm


def plotSunEarth(sunEarth, outfileName, solverType, N, finalTime):
    dt = finalTime/N
    fig, ax = plt.subplots()
    ax.plot(sunEarth.x[:-1], sunEarth.y[:-1])
    ax.set_title(solverType + '\n T %d, $\Delta$ t %.2g' %(finalTime, dt))
    ax.set_xlabel('x [Au]')
    ax.set_ylabel('y [Au]')
    #ax.set_xlim(-2., 2.)
    #ax.set_ylim(-2., 2.)
    fig.savefig('results/' + outfileName + solverType + '.png') 
    #plt.show()
    plt.close()
    

def plotSunEarthTimes(sunEarth, outfileName, solverType, N, finalTime):
    dt = finalTime/N
    fig, ax = plt.subplots(2, sharex=True)
    fig.hold('on')
    ax[0].plot(sunEarth.time[:-1], sunEarth.x[:-1])
    ax[1].plot(sunEarth.time[:-1], sunEarth.y[:-1])
    fig.title(outfileName + ' ' + solverType + '\n N %d, dt %.2g' %(N, dt))
    ax[0].set_xlabel('t [Au]')
    ax[0].set_ylabel('x [Au]')
    ax[1].set_ylabel('y [Au]')
    #ax.set_xlim(-2., 2.)
    fig.ylim(-1.1, 1.1)
    #fig.savefig('results/' + outfileName + 'times.png') 
    #plt.show()
    plt.close()

#%% 2, run 

def sunEarthTerminalVelocity():
    """

    """
    sunEarth      = OrderedDict()
    
    outfileName = 'sunEarthTerminalVelocity'
    solverType = 'VelocityVerlet'
    
    finalTime = 10.**5
    #Ns = [10**i for i in xrange(3,8)]
    dt = 0.01
    N = finalTime/dt
    numberOfSimulations = 3
    epsilon = 1.0E-4
    start = 2.*np.sqrt(2)-epsilon # start value Vy, when Vy = start*2pi/128
    initialVelocities = [(start*np.pi)*(1+epsilon)**i for i in xrange(3)]
    
    # Plots. Positions vs time
    

    # Times
    fig, ax = plt.subplots()
    fig.hold('on')
    ax.set_title(solverType + ' radial distance' + ' $\Delta t$ %.2g' %dt)
    ax.set_xlabel("time")
    ax.set_ylabel('radial distance')
    legends = []
    
    for initialVelocityY in initialVelocities:
        outfileName2 = outfileName + 'initialVelocityY%spi' %str(initialVelocityY/np.pi).replace(".", "")
        print outfileName2 
        runCpp(outfileName2, finalTime, N, solverType, initialVelocityY)
        sunEarth['initialVy%fpi' %(initialVelocityY/np.pi)] = pd.read_table("results/" + outfileName2 + ".csv", 
            			            delimiter=',')
        plotSunEarth(sunEarth['initialVy%fpi' %(initialVelocityY/np.pi)], outfileName2, solverType, N, finalTime)
        ax.plot(sunEarth['initialVy%fpi' %(initialVelocityY/np.pi)].time[:-1], sunEarth['initialVy%fpi' %(initialVelocityY/np.pi)].r[:-1])
        legends.append('initialVy %fpi' %(initialVelocityY/np.pi))
    
    ax.legend(legends, bbox_to_anchor=(1.1, 0.7), ncol=1)
    #ax.set_xlim(-2., 2.)
    #ax.set_ylim(0.9, 1.1)
    fig.savefig('results/'+ outfileName + 'radialDistance' + solverType + '.png')
    plt.close()
    
    return sunEarth

#%% 2, run     
sunearth, supNormValues, supNormAngularMomentum = sunEarth()
#sunearthTerminalVelocity = sunEarthTerminalVelocity()
#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict
from matplotlib.pyplot import cm
from matplotlib import colors as mcolors
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


#%% clean results and movie directory
call(["./Allclean"])       

#%% Run 

def runCpp(outfileName, finalTime, N, solverType, initialVy, beta, scenario):
    """
    Compiles and runs cpp program from the command line and
    makes sure the mesh size is not too big.
    """
    N = str(N)
    finalTime = str(finalTime)
    initialVy = str(initialVy)
    beta = str(beta)
    call(["./AllrunVectorized", outfileName, finalTime, N, solverType, initialVy, beta, scenario])
    
#%% 2, run 

def sunEarth():
    """

    """
    scenario = "twoBody"
    initialVy = 2*np.pi
    sunEarth      = OrderedDict()
    supNormValues = []#OrderedDict()
    supNormAngularMomentum = []
    timesUsed = []
    beta = 3.0
    outfileName = 'sunEarth'
    solverType = 'VelocityVerlet'
    
    finalTimes = [10**i for i in xrange(0,3)]
    #Ns = [10**i for i in xrange(3,8)]   
    dts = [10.**(-i) for i in xrange(1, 4)]
    
    Ns = np.asarray(finalTimes)/np.asarray(dts)
    
    
    # Plots. Positions vs time
    fig, ax = plt.subplots(2, sharex=True)
    fig.hold('on')
    ax[0].set_title(outfileName + ' ' + solverType)
    ax[1].set_xlabel('t [Au]')
    ax[0].set_ylabel('x [Au]')
    ax[1].set_ylabel('y [Au]')
    #plt.ylim(-1.1, 1.1)    
    legends = []
    
    # Energy plots
    fig2, ax2 = plt.subplots()
    fig2.hold('on')
    ax2.set_title(outfileName + ' ' + solverType + 'Energy')
    ax2.set_xlabel('t [Au]')
    ax2.set_ylabel('E ')
    #ax2.set_ylim(0.7, 1.05)    
    

    
    # Angular momentum plots
    fig3, ax3 = plt.subplots()
    fig3.hold('on')
    ax3.set_title(outfileName + ' ' + solverType + 'Angular momentum')
    ax3.set_xlabel('t [Au]')
    ax3.set_ylabel('Angular momentum ')
    #ax3.set_ylim(0.99, 1.05)    
    
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
            print 'N = %d, final time = %.2g ' %(N, finalTime)  
            outfileName2 = outfileName + 'finalTime%s' %str(finalTime).replace(".", "") + 'N%s' %str(int(round(N)))#.replace(".", "")
            runCpp(outfileName2, finalTime, N, solverType, initialVy, beta, scenario)
            sunEarth['FinalTime %f' %finalTime]['N %f' %N] = pd.read_table("results/" + outfileName2 + "Earth.csv", 
            			            delimiter=',')
            plotSunEarth(sunEarth['FinalTime %f' %finalTime]['N %f' %N], outfileName2, solverType, N, finalTime)
            if finalTime == finalTimes[1]:
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
    #fig.ylim(-1.1, 1.1)
    #fig.savefig('results/' + outfileName + 'times.png') 
    #plt.show()
    plt.close()

#%% 2, run 

def sunEarthTerminalVelocity():
    """

    """
    scenario = "twoBody"
    sunEarth      = OrderedDict()
    
    outfileName = 'sunEarthTerminalVelocity'
    solverType = 'VelocityVerlet'
    beta = 2.0
    finalTime = 10.**3
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
        runCpp(outfileName2, finalTime, N, solverType, initialVelocityY, beta, scenario)
        sunEarth['initialVy%fpi' %(initialVelocityY/np.pi)] = pd.read_table("results/" + outfileName2 + "Earth.csv", 
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

#%% 
def sunEarthAlternativeGravitationalForce():
    """

    """
    sunEarth      = OrderedDict()
    
    outfileName = 'sunEarth'
    solverType = 'AlternativeForce'
    betaValues = [3.0,  3.5, 3.9, 4.0]
    finalTime = 10.**4
    #Ns = [10**i for i in xrange(3,8)]
    dt = 0.01
    N = finalTime/dt
    numberOfSimulations = 5
    epsilon = np.pi/8
    start = 2.*np.pi # start value Vy, when Vy = start*2pi/128
    initialVelocities = [start+(epsilon)*i for i in xrange(numberOfSimulations)]
    
    # Plots. Positions vs time
    
    for beta in betaValues:
        sunEarth['beta %f' %beta] = {}
        fig, ax = plt.subplots()
        fig.hold('on')
        #ax.set_title(solverType + ' radial distance' + ' $\Delta t$ %.2g' %dt + "$\beta $" + str(beta))
        ax.set_xlabel("time")
        ax.set_ylabel('radial distance')
        legends = []
        
        if abs(beta - 4.0) < 1E-12:
            initialVelocities = [2*np.pi - 1*np.pi/1024 + i*np.pi/1024 for i in xrange(3)] 
        for initialVelocityY in initialVelocities:
        
            outfileName2 = outfileName + 'initialVelocityY%spibeta%s' %(str(initialVelocityY/np.pi).replace(".", ""), str(beta).replace(".", ""))
            print outfileName2 
            runCpp(outfileName2, finalTime, N, solverType, initialVelocityY, beta, scenario)
            sunEarth['beta %f' %beta]['initialVy%fpi' %(initialVelocityY/np.pi)] = pd.read_table("results/" + outfileName2 + ".csv", 
                			            delimiter=',')
            plotSunEarth(sunEarth['beta %f' %beta]['initialVy%fpi' %(initialVelocityY/np.pi)], outfileName2, solverType, N, finalTime)
            ax.plot(sunEarth['beta %f' %beta]['initialVy%fpi' %(initialVelocityY/np.pi)].time[:-1], sunEarth['beta %f' %beta]['initialVy%fpi' %(initialVelocityY/np.pi)].r[:-1])
            legends.append('initialVy %fpi' %(initialVelocityY/np.pi))
        
        ax.legend(legends, bbox_to_anchor=(1.1, 0.7), ncol=1)
        #ax.set_xlim(-2., 2.)
        #ax.set_ylim(0.9, 1.1)
        fig.savefig('results/'+ outfileName + 'radialDistance' + solverType + 'beta' + str(beta).replace(".", "")+'.png')
        plt.close()
    
    return sunEarth

#%% Multibody, sun stationar     

def multiBodyStationarySun(threePlanets, threeBodiesMovingSun,scenario,movie):
    """

    """
    if not os.path.isdir('movie'):
        os.mkdir('movie')
    if threePlanets:
        planets = ['Earth', 'Jupiter']
    if threeBodiesMovingSun:
        planets = ['Sun','Earth', 'Jupiter']
        
    elif scenario == 'solarSystem':
        planets = ['Sun', 'Earth', 'Jupiter', 'Mars', 'Venus', 'Saturn', 'Mercury', 'Uranus', 'Neptune', 'Pluto']
    elif scenario == 'mercury':
        planets = ['Mercury']

    initialVy = 2*np.pi
    beta = 3.0
    outfileName = scenario
    solverType = 'VelocityVerlet'
    
    finalTimes = [10**i for i in xrange(2,3)]
    #Ns = [10**i for i in xrange(3,8)]   
    dts = [10.**(-i) for i in xrange(7, 8)]
    
    Ns = np.asarray(finalTimes)/np.asarray(dts)
    
    multiBodies      = OrderedDict()
    for finalTime in finalTimes:
        multiBodies['FinalTime %f' %finalTime] = {}
        for dt in dts:
            multiBodies['FinalTime %f' %finalTime]['dt %f' %dt] = {}
            N = finalTime/dt
            print 'N = %d, final time = %.2g ' %(N, finalTime)  
            outfileName2 = outfileName + 'finalTime%s' %str(finalTime).replace(".", "") + 'N%s' %str(int(round(N)))#.replace(".", "")
            runCpp(outfileName2, finalTime, N, solverType, initialVy, beta, scenario)
            fig, ax = plt.subplots()
            legends = []
            plt.hold('on')
            
            for planet in planets:
                multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet] = pd.read_table("results/" + outfileName2 + planet + ".csv", 
            			            delimiter=',')
                
                ax.plot(multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].x[:-1], multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].y[:-1])
                legends.append(planet)
                
            plt.legend(legends)
            dt = finalTime/N
            ax.set_title(solverType + '\n T %d, $\Delta$ t %.2g' %(finalTime, dt))
            ax.set_xlabel('x [Au]')
            ax.set_ylabel('y [Au]')
            plt.axis('equal')
            #ax.set_xlim(-2., 2.)
            #ax.set_ylim(-2., 2.)
            fig.savefig('results/' + outfileName + solverType + 'T' + str(finalTime).replace(".", "") + 'dt' + str(dt).replace(".", "") + '.png') 
            plt.close()
            
            # Movie
            if movie:
                if finalTime == 100 and dt == 0.0001:
                    saveInterval = 200
                    numberOfObservations = len(multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].time[:-1]) 
                    numberOfPlots = int(round(numberOfObservations/saveInterval))
                    fig2, ax2 = plt.subplots()
                    plt.hold('on')
                    #ax2.set_title(solverType + '\n T %d, $\Delta$ t %.2g Time between frames %.1g' %(finalTime, dt, dt*saveInterval))
                    ax2.set_xlabel('x [Au]')
                    ax2.set_ylabel('y [Au]')
                    plt.axis('equal')
                    ax2.set_xlim(-2., 2.)
                    ax2.set_ylim(-2., 2.)
    
                    colors = ['b', 'g', 'y', 'r', 'm' ,'b' , 'y', 'm', 'r']#, 'coral', 'cornsilk']
                    #colors = ['black', 'red', 'green']#, 'palegreen', 'blue', 'mediumturquoise', 'deeppink', 'purple', 'yellow']
                    #color=iter(cm.rainbow(np.linspace(0,1,9)))
                    """for i in range(n):
                       c=next(color)
                       ax1.plot(x, y,c=c)"""
                    #colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive']
                    for counter in xrange(numberOfPlots):
                        for planet,color in zip(planets, colors):
                            ax2.plot(multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].x[counter*saveInterval], multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].y[counter*saveInterval], 'o%s'%color, ms=5)#  
                            ax2.set_title(solverType + '\n T %d, $\Delta$ t %.2g Time between frames %.1g \n t %.2g' %(finalTime, dt, dt*saveInterval, multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].time[counter*saveInterval]))
                            
    
                        ax2.legend(legends, bbox_to_anchor=(1.05, 0.7), ncol=1)
                        fig2.savefig('movie/tmp_%04d' %counter + outfileName + solverType + 'T' + str(finalTime).replace(".", "") + 'dt' + str(dt).replace(".", "") + '.png') 
                    plt.close()
                    
                    # Make video file
                    fps = 4  # frames per second
                    codec2ext = dict(libx264='mp4')  # video formats
                    filespec = 'movie/tmp_%04d' + outfileName + solverType + 'T' + str(finalTime).replace(".", "") + 'dt' + str(dt).replace(".", "") + '.png'
                    movie_program = 'ffmpeg'  # or 'avconv'
                    for codec in codec2ext:
                        ext = codec2ext[codec]
                        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s ' \
                              '-vcodec %(codec)s movie3_6.%(ext)s' % vars()
                        os.system(cmd)
                        
            # Perillion precession
            if scenario == 'mercury':
                rPerrilion = np.min(multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].r.values)
                indexRPerrilion = np.argmin(multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].r.values)
                xPerrilion = multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].x.values[indexRPerrilion]
                yPerrilion = multiBodies['FinalTime %f' %finalTime]['dt %f' %dt]['Planet %s' %planet].y.values[indexRPerrilion]
                theta = np.arctan(yPerrilion/xPerrilion)
                precession = theta*206265.
                print 'Perillion precessiun %.4f' %precession
                print 'theta %.4f' %theta
                
    return multiBodies

#%% 2, run     
#sunearth, supNormValues, supNormAngularMomentum = sunEarth()
#sunearthTerminalVelocity = sunEarthTerminalVelocity()
#sunEarthAlternativeGravitationalForce = sunEarthAlternativeGravitationalForce()
#miltiBodies = multiBodyStationarySun(threePlanets = True, threeBodiesMovingSun=False, scenario = "threeBodiesJupiterTimes1000")
#miltiBodies = multiBodyStationarySun(threePlanets = False,threeBodiesMovingSun=True, scenario = "threeBodiesJupiterMassTimes1000MovingSun")
#miltiBodies = multiBodyStationarySun(threePlanets = False,threeBodiesMovingSun=False, scenario = "solarSystem")
miltiBodies = multiBodyStationarySun(threePlanets = False,threeBodiesMovingSun=False, scenario = 'mercury', movie=False)


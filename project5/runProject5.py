#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy.misc import comb
from mpl_toolkits.mplot3d import Axes3D  # necessary for 3D plotting
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#from scitools.std import *
import scitools.std as st
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL) 

#%% Run 
class Project5:
    
    def __init__(self):
        return None
        
    def runCpp(self, outfileName, dt, dx, theta, T, scenario, threadNumber):
        """
        Compiles and runs cpp program from the command line and
        makes sure the mesh size is not too big.
        """
        dt = str(dt)
        dx = str(dx)
        theta = str(theta)
        T = str(T)
        threadNumber = str(threadNumber)
        call(["./AllrunVectorized", outfileName, dt, dx, theta, T, scenario, threadNumber])
        
    def project5c(self):
        dxValues = [0.1]#, 0.01]
        safetyFacors = [0.96, 1.04]
        movieCounter = 0
        dimension = "1D"
        threadNumber = 8
        for dx in dxValues:
            for safetyFactor in safetyFacors:
                movieCounter += 1
                dt = dx**2/2.0*(1/safetyFactor)
                alpha = dt/dx**2
                theta = 0.5
                T = 0.5
                nT = T/dt+1
                outfileName ='out5C'
                outfileName2 = os.getcwd() + '/results/' + outfileName
                saveEveryNSolution = 10
                
                self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)
                
                # Read data
                data = OrderedDict()
                scenerios = ['BackwardEuler', 'ForwardEuler', 'CrancNicholson', 'Analytical']
                x = np.zeros(int(round(1./dx+1)))
                t = np.zeros(int(round(T/dt +1)))
                for scenario in scenerios:    
                    data[scenario] = pd.read_csv(outfileName2 + scenario + 'SolutionMatrixU.txt', delim_whitespace=True, header=None)
                x = data[scenario].ix[1:,0]
                t = data[scenario].ix[0,1:]
                
                #results4b.columns = ["acceptedMoves","mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                
                #Plots
                counter =1
                figureNumber = 1
                fig2,ax2 = plt.subplots()
                ax2.hold(True)
                ax2Title = []
                
                for times in t:
                    fig,ax = plt.subplots()
                    ax.set_xlim([0, 1])
                    ax.set_ylim([0, 1])
                    ax.hold(True)
                    if counter%saveEveryNSolution == 0 or counter == 1:
                        for scenario in scenerios:
                             ax.plot(x,data[scenario].ix[1:,counter])
                        figureNumber += 1        
                        ax.set_xlabel(r'$x$')
                        ax.set_ylabel(r"$u(x,t)$")
                        ax.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %g $ ' %alpha + '\n' + 't = %.5g' %times+ '\n')
                        ax.legend(scenerios, loc=2)
                        fig.tight_layout()
                        fig.savefig('movie/%1dtmp_%04d' %(movieCounter, figureNumber) + outfileName + '.png') 
                        plt.close()
                    counter += 1
        
                    if counter == 2 or counter == len(t)/6 or counter == len(t):
                        for scenario in scenerios:
                            if scenario == 'ForwardEuler':
                                markerType = 'r-+'
                            elif scenario == 'BackwardEuler':
                                markerType = 'b-o'
                            elif scenario == 'CrancNicholson':
                                markerType = 'g-^'
                            else:
                                markerType = 'y-'
                            ax2.plot(x,data[scenario].ix[1:,counter], markerType)
                        ax2Title.append(round(times,4))
        
                ax2.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %g $ ' %alpha + '\n' + 't = %s' %ax2Title + '\n')
                ax2.legend(scenerios, loc=2)
                ax2.set_xlabel(r'$x$')
                ax2.set_ylabel(r"$u(x,t)$")
                fig2.tight_layout()
                fig2.savefig(outfileName2 + 'Number' + str(movieCounter) + '.png') 
                plt.close()
        
                
#                # Make video file
#                fps = 4  # frames per second
#                codec2ext = dict(libx264='mp4')  # video formats
#                filespec = 'movie/%1dtmp_%04d' %movieCounter + outfileName + '.png'
#                movie_program = 'ffmpeg'  # or 'avconv'
#                for codec in codec2ext:
#                    ext = codec2ext[codec]
#                    cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s ' \
#                          '-vcodec %(codec)s movie3_6.%(ext)s' % vars()
#                    os.system(cmd)
        return data
       
        
    def project5d(self):
        dxValues = [0.1, 0.05, 0.025, 0.0125]#, 0.00625]#, 0.025]#, 0.01]#, 0.01]
        safetyFactor = 1.05
        dt = dxValues[-1]**2/2.0*(1/safetyFactor)
        T = 0.25
        theta = 0.5
        dimension = "1D"
        threadNumber = 8
        outfileName ='out5NormDx'
        outfileName2 = os.getcwd() + '/results/' + outfileName
        
        counter = 1
        data = OrderedDict()
        for dx in dxValues:
            alpha = dt/dx**2
            self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)            
            # Read data
            data[counter] = pd.read_csv(outfileName2 + '.txt', delim_whitespace=True, header=None) 
            counter += 1
            
        #Plots
        numberOfScenarios = 3
        numerOfDx = len(dxValues)
        scenarioNr = range(1,numerOfDx+1)
        scenerios = ['BackwardEuler', 'ForwardEuler', 'CrancNicholson']#, 'Analytical']
        errorFe = [data[i].ix[0,2] for i in scenarioNr]
        errorBe = [data[i].ix[1,2] for i in scenarioNr]
        errorCn = [data[i].ix[2,2] for i in scenarioNr]
        errors = [errorFe, errorBe, errorCn]
        
        fig2,ax2 = plt.subplots()
        logDxValues = []
        for i in xrange(numerOfDx):
            logDxValues.append(np.log2(dxValues[i])) 

        for error in errors:
            ax2.plot(logDxValues, error, '-o')
        ax2.set_xlabel(r'$\log_2 \Delta x$')
        ax2.set_ylabel(r"$\log_2 |Error|_{\infty}$")
        ax2.set_title('Sup norm x')
        ax2.legend(scenerios, loc=2) 
        fig2.tight_layout()
        fig2.savefig(outfileName2+ '.png')
        plt.close()
        
        return data

    def project5f(self):
            dxValues = [0.1]#, 0.01]
            safetyFacors = [1.04]#[0.96, 1.04]
            movieCounter = 0
            dimension = "2DExplicit"
            threadNumber = 8
            for dx in dxValues:
                nx = int(round(1./dx + 1))
                x = np.linspace(0,1,nx)
                y = x
                X, Y = np.meshgrid(x, y)
                for safetyFactor in safetyFacors:
#                    fig = plt.figure()
#                    ax = fig.gca(projection='3d')
                    
                    #cv = u_box.grid.coorv  # vectorized mesh coordinates
                    movieCounter += 1
                    dt = dx**2/4.0*(1/safetyFactor)
                    alpha = dt/dx**2
                    theta = 0.5
                    T = .3
                    nT = int(round(T/dt+1))
                    outfileName ='out5f'
                    outfileName2 = os.getcwd() + '/results/' + outfileName
                    saveEveryNSolution = 10
                    
                    self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)
                    
                    # Read data
                    #data = OrderedDict()
                    scenerios = ['Explicit', 'Analytical']
                    #x = np.zeros(int(round(1./dx+1)))
                    #t = np.zeros(int(round(T/dt +1)))
                    xv = x.reshape((x.size,1)) #x[:, np.newaxis]          # for vectorized function evaluations
                    yv = y.reshape((1,y.size))#y[np.newaxis, :]
                        #data[scenario] = OrderedDict()
                    fileCounter = 1
                    for fileCounter in xrange(1, nT):
                        data = pd.read_csv(outfileName2 + 'SolutionMatrixUTime%d.txt' %fileCounter, delim_whitespace=True, header=None)
                        
                        #surf = ax.plot_surface(X, Y, data.values, cmap=cm.coolwarm,
                        #rstride=1, cstride=1, linewidth=0, antialiased = False)
                        
                        st.surfc(xv, yv, data, title='Explicit Time = %.4f' %((fileCounter-1)*dt) , zlim=[-0.1, 1.1],
                              colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                              shading='flat')  #
                        st.savefig('movie/tmp_%04d' %fileCounter + outfileName + '.png') 
                        
                        data2 = pd.read_csv(outfileName2 + 'AnalyticalSolutionMatrixU2D%d.txt' %fileCounter, delim_whitespace=True, header=None)
                        st.surfc(xv, yv, data2, title='Analytical Time = %.4f' % ((fileCounter-1)*dt), zlim=[-0.1, 1.1],
                              colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                              shading='flat')  #
                        st.savefig('movie/tmpAnalytic_%04d' %fileCounter + outfileName + '.png')
                        
                        st.surfc(xv, yv, data-data2, title='Numerical-Analytical time = %.4f' % ((fileCounter-1)*dt), zlim=[-.1, 1.1],
                              colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                              shading='flat')  #
                        st.savefig('movie/tmpError_%04d' %fileCounter + outfileName + '.png')
                        
                        
            
            data = pd.read_csv(outfileName2 + 'Timing.txt' , delimiter=',')
               
            width = 0.35 
            fig, ax = plt.subplots()
            rects1 = ax.bar(np.arange(3), [np.asscalar(data.numerical.values), np.asscalar(data.analytic.values), np.asscalar(data.analytic.values)/np.asscalar(data.numerical.values)], width)#, color='r')
            fontSizes = 20
            ax.set_title('Timing Numerical and analytical', fontsize = fontSizes )
            ax.set_xticks(np.arange(3) + width / 10)
            ax.set_xticklabels(('Numerical', 'Analytical', r'$\frac{Analytical}{Numerical}$'))
            #plt.show()
            fig.tight_layout()
            plt.savefig( outfileName2 + 'Timing.png')
            return data
        
    def project5g(self):
        dx = 0.1#, 0.01]
        safetyFactor = 1.04#[0.96, 1.04]
        dimension = "2DExplicit"
        threadNumbers = [1,8]
        dt = dx**2/4.0*(1/safetyFactor)
        alpha = dt/dx**2
        theta = 0.5
        T = .1
        nT = int(round(T/dt+1))
        
        outfileName ='out5g'
        outfileName2 = os.getcwd() + '/results/' + outfileName
        
        data = OrderedDict()
        
        for threadNumber in threadNumbers: 
            self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)
            data[threadNumber] = pd.read_csv(outfileName2 + 'Timing.txt' , delimiter=',')
        
        width = 0.35 
        fig, ax = plt.subplots()
        rects1 = ax.bar(np.arange(2), [np.asscalar(data[1].analytic.values), np.asscalar(data[8].analytic.values)], width)#, color='r')
        fontSizes = 20
        ax.set_title('Timing parallelization analytical', fontsize = fontSizes )
        ax.set_xticks(np.arange(2) + width / 10)
        ax.set_xticklabels(('Unparallelized', 'Parallelized'))
        #plt.show()
        fig.tight_layout()
        plt.savefig( outfileName2 + 'TimingAnalytical.png')
        return data
#%%
if __name__ == "__main__":
    # Clean results and movie directory
    call(["./Allclean"]) 
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
        
    if not os.path.isdir('movie'):
        os.mkdir('movie')  
    
    parser = argparse.ArgumentParser(description="starts a c++ program simulating heat equation, generatsplots. \n Alternatives: 5c")
    parser.add_argument("task", type=str, default='5c', help="press 5c")
    args = parser.parse_args()
    
    if args.task  == '5c':
       project5c = Project5()
       data = project5c.project5c()

    if args.task  == '5d':
       project5d = Project5()
       data=project5d.project5d()

    if args.task  == '5f':
       project5f = Project5()
       data = project5f.project5f()
    
    if args.task  == '5g':
       project5g = Project5()
       data = project5g.project5g()
       
#%%



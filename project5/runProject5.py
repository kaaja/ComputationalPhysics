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

#%% Run 
class Project5:
    
    def __init__(self):
        return None
        
    def runCpp(self, outfileName, dt, dx, theta, T):
        """
        Compiles and runs cpp program from the command line and
        makes sure the mesh size is not too big.
        """
        dt = str(dt)
        dx = str(dx)
        theta = str(theta)
        T = str(T)
        call(["./AllrunVectorized", outfileName, dt, dx, theta, T])
        
    def project5c(self):
        dxValues = [0.1]#, 0.01]
        safetyFacors = [0.96, 1.04]
        movieCounter = 0
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
                
                self.runCpp(outfileName2, dt, dx, theta, T)
                
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
                fig2.savefig(outfileName + 'Number' + str(movieCounter) + '.png') 
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
        
    def project5d(self):
        dxValues = [0.1, 0.05]#, 0.01]#, 0.01]
        safetyFactor = 1.05
        dt = dxValues[-1]**2/2.0*(1/safetyFactor)
        theta = 0.5
        T = 0.5
        
        counter = 1
        data = OrderedDict()
        fig2,ax2 = plt.subplots()
        ax2.hold(True)
        for dx in dxValues:
            alpha = dt/dx**2
                
            nT = T/dt+1
            outfileName ='out5d' + 'dx%s' %str(counter)
            outfileName2 = os.getcwd() + '/results/' + outfileName
            
            self.runCpp(outfileName2, dt, dx, theta, T)
            
            # Read data
           
            scenerios = ['BackwardEuler', 'ForwardEuler', 'CrancNicholson', 'Analytical']
                
            data = pd.read_csv(outfileName2 + '.txt', delim_whitespace=True, header=None) 
            logdx.append(data.ix[:,1])
            logdt.append(data.ix[:,1])
            logError.append(data.ix[:,2])
            
            
            
            counter += 1
            
            
            #results4b.columns = ["acceptedMoves","mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
        
        #Plots
        ax2Title = []
        
        ax2.plot([data['1'].ix[1,1],data['2'].ix[1,1]], [data['1'].ix[0,2],data['2'].ix[0,2]])
#        figureNumber += 1        
#        ax.set_xlabel(r'$x$')
#        ax.set_ylabel(r"$u(x,t)$")
#        ax.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %g $ ' %alpha + '\n' + 't = %.5g' %times+ '\n')
#        ax.legend(scenerios, loc=2)
#        fig.tight_layout()
        fig2.savefig(outfileName2+ '.png') 
        plt.close()
    
#            if counter == 2 or counter == len(t)/6 or counter == len(t):
#                for scenario in scenerios:
#                    if scenario == 'ForwardEuler':
#                        markerType = 'r-+'
#                    elif scenario == 'BackwardEuler':
#                        markerType = 'b-o'
#                    elif scenario == 'CrancNicholson':
#                        markerType = 'g-^'
#                    else:
#                        markerType = 'y-'
#                    ax2.plot(x,data[scenario].ix[1:,counter], markerType)
#                ax2Title.append(round(times,4))
#            counter += 1
#    
#        ax2.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %g $ ' %alpha + '\n' + 't = %s' %ax2Title + '\n')
#        ax2.legend(scenerios, loc=2)
#        ax2.set_xlabel(r'$x$')
#        ax2.set_ylabel(r"$u(x,t)$")
#        fig2.tight_layout()
#        fig2.savefig(outfileName + 'Number' + str(movieCounter) + '.png') 
#        plt.close()
    
            
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
       project5c.project5c()

    if args.task  == '5d':
       project5d = Project5()
       data=project5d.project5d()

       
#%%



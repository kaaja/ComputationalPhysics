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
        dxValues = [0.1, 0.01]
        safetyFacors = [0.96, 1.04]
        movieCounter = 0
        dimension = "1D"
        threadNumber = 8
        showMovie = False
        ForwardEulerInFigures = False
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
                scenariosWithoutForwardEuler = ['BackwardEuler', 'CrancNicholson', 'Analytical']
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
                #ax2.hold(True)
                ax2Title = []
                
                for times in t:
                    if showMovie:
                        fig,ax = plt.subplots()
                        ax.set_xlim([0, 1])
                        ax.set_ylim([0, 1])
                        #ax.hold(True)
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
        
                    if ForwardEulerInFigures:
                            scenariosLoop = scenarios
                    else:
                            scenariosLoop = scenariosWithoutForwardEuler
        
                    if  counter == len(t)/30 or counter == len(t)/3 :
                        for scenario in scenariosLoop:
                            if scenario == 'ForwardEuler':
                                markerType = 'r-'
                            elif scenario == 'BackwardEuler':
                                markerType = 'b-'
                            elif scenario == 'CrancNicholson':
                                markerType = 'g-'
                            else:
                                markerType = 'y-'
                            ax2.plot(x,data[scenario].ix[1:,counter], markerType)
                        ax2Title.append(round(times,4))

        
                ax2.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %g $ ' %alpha + '\n' + 't = %s' %ax2Title + '\n')
                ax2.legend(scenariosLoop, loc=2)
                ax2.set_xlabel(r'$x$')
                ax2.set_ylabel(r"$u(x,t)$")
                fig2.tight_layout()
                fig2.savefig(outfileName2 + 'Number' + str(movieCounter) + '.png') 
                plt.close()
                #plt.show()
                
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
        #dx = 0.1#, 0.05, 0.025, 0.0125]#, 0.00625]#, 0.025]#, 0.01]#, 0.01]
        
        #dtStart = .005
        #dtStart = dtStart*(1./2)**7
        #dtStart = dtStart*(1./2)**5
        #numberOfDts = 2
        #dtValues = [dtStart*(1./2)**i for i in xrange(numberOfDts)]
        
        safetyFacors = [0.96, 1.04]
        dxValues = [0.1, 0.01]
        times1 = [0.048 ,0.1587, 0.2, 0.25]
        times2 = [0.048,0.1587]
        #dt = dxValues[-1]**2/2.0*(1/safetyFactor)
        theta = 0.5
        dimension = "1D"
        threadNumber = 8
        outfileName ='out5NormDx'
        outfileName2 = os.getcwd() + '/results/' + outfileName
        
        
        outfileNameCounter = 1
        for dx in dxValues:
            for safetyFactor in safetyFacors:
                if dx == 0.01 and safetyFactor == 0.96:
                    endTimeValues = times2
                else:
                    endTimeValues = times1
                counter = 1
                data = OrderedDict()
                for T in endTimeValues:
                    #dx = dt
                    dt = dx**2/2.0*(1/safetyFactor)
                    alpha = dt/dx**2
                    self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)            
                    # Read data
                    data[counter] = pd.read_csv(outfileName2 + '.txt', delim_whitespace=True, header=None) 
                    counter += 1
                    
                #Plots
                numberOfScenarios = 3
                numerOfEndTimes = len(endTimeValues)
        
                numerOfDt = len(endTimeValues)
                scenarioNr = range(1,numerOfEndTimes+1)
                scenerios = ['Forward Euler', 'Backward Euler','Crank-Nicolson']#, 'Analytical']
                errorFe = [data[i].ix[0,2] for i in scenarioNr]
                errorBe = [data[i].ix[1,2] for i in scenarioNr]
                errorCn = [data[i].ix[2,2] for i in scenarioNr]
                errors = [np.asarray(errorFe), np.asarray(errorBe), np.asarray(errorCn)]
                
                fig2,ax2 = plt.subplots()
                #logDxValues = []
        
                ind = np.arange(numerOfEndTimes)  # the x locations for the groups
                width = 0.25       # the width of the bars
                rects1 = ax2.bar(ind, errors[0], width, color='r')#, yerr=men_std)
                rects2 = ax2.bar(ind+width, errors[1], width, color='b')#, yerr=men_std)
                rects3 = ax2.bar(ind+2*width, errors[2], width, color='g')#, yerr=men_std)
                ax2.set_xlabel(r'$Time$')
                ax2.set_ylabel(r"$L_2$ norm")
                if dx == 0.1:
                    ax2.set_ylim(0,0.0150)
                ax2.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %.2g $ ' %alpha + '\n')
                ax2.set_xticks(ind + width / 2)
                ax2.set_xticklabels(endTimeValues)
                ax2.legend(scenerios, loc=2) 
                #ax2.legend(scenerios,bbox_to_anchor=(1.04,1), borderaxespad=0)              
                #fig2.tight_layout(rect=[0,0,0.75,1])
                fig2.tight_layout()
                fig2.savefig(outfileName2+ str(outfileNameCounter) + '.png')
                #fig2.savefig(outfileName2+ str(outfileNameCounter) + '.png', bbox_inches="tight")
                #plt.show()
                plt.close()
                outfileNameCounter+= 1
        
        return data
    
    def project5d2D(self):
        #dx = 0.1#, 0.05, 0.025, 0.0125]#, 0.00625]#, 0.025]#, 0.01]#, 0.01]
        
        #dtStart = .005
        #dtStart = dtStart*(1./2)**7
        #dtStart = dtStart*(1./2)**5
        #numberOfDts = 2
        #dtValues = [dtStart*(1./2)**i for i in xrange(numberOfDts)]
        
        safetyFacors = [0.8, 1.04]
        dxValues = [0.1]#, 0.01]
        times1 = [0.0125, 0.0625]#, 0.1167]
        times2 = [0.0125, 0.0625]#, 0.1167]
        #dt = dxValues[-1]**2/2.0*(1/safetyFactor)
        theta = 0.5
        dimension = "2D"
        threadNumber = 8
        outfileName ='out5Norm2D'
        outfileName2 = os.getcwd() + '/results/' + outfileName
        
        
        outfileNameCounter = 1
        for dx in dxValues:
            for safetyFactor in safetyFacors:
                if dx == 0.01 and safetyFactor == 0.96:
                    endTimeValues = times2
                else:
                    endTimeValues = times1
                counter = 1
                data = OrderedDict()
                for T in endTimeValues:
                    #dx = dt
                    dt = dx**2/4.0*(1/safetyFactor)
                    alpha = dt/dx**2
                    self.runCpp(outfileName2, dt, dx, theta, T, dimension, threadNumber)            
                    # Read data
                    data[counter] = pd.read_csv(outfileName2 + '.txt', delim_whitespace=True, header=None) 
                    counter += 1
                    
                #Plots
                numberOfScenarios = 3
                numerOfEndTimes = len(endTimeValues)
        
                numerOfDt = len(endTimeValues)
                scenarioNr = range(1,numerOfEndTimes+1)
                scenerios = ['Forward Euler', 'Backward Euler']#,'Crank-Nicolson']#, 'Analytical']
                errorFe = [data[i].ix[0,2] for i in scenarioNr]
                errorBe = [data[i].ix[1,2] for i in scenarioNr]
                #errorCn = [data[i].ix[2,2] for i in scenarioNr]
                errors = [np.asarray(errorFe), np.asarray(errorBe)]#, np.asarray(errorCn)]
                
                fig2,ax2 = plt.subplots()
                #logDxValues = []
        
                ind = np.arange(numerOfEndTimes)  # the x locations for the groups
                width = 0.25       # the width of the bars
                rects1 = ax2.bar(ind, errors[0], width, color='r')#, yerr=men_std)
                rects2 = ax2.bar(ind+width, errors[1], width, color='b')#, yerr=men_std)
                #rects3 = ax2.bar(ind+2*width, errors[2], width, color='g')#, yerr=men_std)
                ax2.set_xlabel(r'$Time$')
                ax2.set_ylabel(r"$L_2$ norm")
                '''if dx == 0.1:
                    ax2.set_ylim(0,0.0150)'''
                ax2.set_title(r'$\Delta x = %.g $' %dx + r' $\frac{\Delta t}{\Delta x^2} = %.2g $ ' %alpha + '\n')
                ax2.set_xticks(ind + width / 2)
                ax2.set_xticklabels(endTimeValues)
                ax2.legend(scenerios, loc=0) 
                #ax2.legend(scenerios,bbox_to_anchor=(1.04,1), borderaxespad=0)              
                #fig2.tight_layout(rect=[0,0,0.75,1])
                fig2.tight_layout()
                fig2.savefig(outfileName2+ str(outfileNameCounter) + '.png')
                #fig2.savefig(outfileName2+ str(outfileNameCounter) + '.png', bbox_inches="tight")
                #plt.show()
                plt.close()
                outfileNameCounter+= 1
        
        return data

    def project5f(self):
            dxValues = [0.1]#, 0.01]
            safetyFacors = [1.04] #[0.8]
            movieCounter = 0
            dimension = "2D"
            threadNumber = 8
            FontSizeUniversal = 22
            showMovie = True
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
                    T = .15
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
                    st.setp(interactive=False)
                    if showMovie:
                        for fileCounter in xrange(1, nT):
                            data = pd.read_csv(outfileName2 + 'ImplicitSolutionMatrixUTime%d.txt' %fileCounter, delim_whitespace=True, header=None)
                            
                            #surf = ax.plot_surface(X, Y, data.values, cmap=cm.coolwarm,
                            #rstride=1, cstride=1, linewidth=0, antialiased = False)
                            
                            st.surfc(xv, yv, data, title='Implicit Time = %.4f' %((fileCounter-1)*dt) , zlim=[-0.1, 1.1],
                                  colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                                  shading='flat', xlabel='x', ylabel='y', zlabel='u')  #
                            st.savefig('movie/tmpImplicit_%04d' %fileCounter + outfileName + '.png')
                            
                            data2 = pd.read_csv(outfileName2 + 'ExplicitSolutionMatrixUTime%d.txt' %fileCounter, delim_whitespace=True, header=None)
                            
                            #surf = ax.plot_surface(X, Y, data.values, cmap=cm.coolwarm,
                            #rstride=1, cstride=1, linewidth=0, antialiased = False)
                            
                            st.surfc(xv, yv, data2, title='Explicit Time = %.4f' %((fileCounter-1)*dt) , zlim=[-0.1, 1.1],
                                  colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                                  shading='flat', xlabel='x', ylabel='y', zlabel='u')  #
                            st.savefig('movie/tmpExplicit_%04d' %fileCounter + outfileName + '.png')  
                            
                            data3 = pd.read_csv(outfileName2 + 'AnalyticalSolutionMatrixU2D%d.txt' %fileCounter, delim_whitespace=True, header=None)
                            st.surfc(xv, yv, data3, title='Analytical Time = %.4f' % ((fileCounter-1)*dt), fontsize = 30,zlim=[-0.1, 1.1],
                                  colorbar=True, colormap=st.hot(), caxis=[-0.1, 1.1],
                                  shading='flat', xlabel='x', ylabel='y', zlabel='u')  #
                            st.savefig('movie/tmpAnalytic_%04d' %fileCounter + outfileName + '.png')
                            
                            st.surfc(xv, yv, (data/data3-1)*100, title='Implicit-Analytical time = %.4f' % ((fileCounter-1)*dt), zlim=[-100.1, 100.1],
                                  colorbar=True, colormap=st.hot(), caxis=[-100.1, 100.1],
                                  shading='flat', xlabel='x', ylabel='y', zlabel='u')  #
                            st.savefig('movie/tmpErrorImplicit_%04d' %fileCounter + outfileName + '.png')
                            
                            st.surfc(xv, yv, (data2/data3-1)*100, title='Explicit-Analytical time = %.4f' % ((fileCounter-1)*dt), zlim=[-100.1, 100.1],
                                  colorbar=True, colormap=st.hot(), caxis=[-100.1, 100.1],
                                  shading='flat', xlabel='x', ylabel='y', zlabel='u')  #
                            st.savefig('movie/tmpErrorExplicit_%04d' %fileCounter + outfileName + '.png')
                            
                        
            
            data = pd.read_csv(outfileName2 + 'Timing.txt' , delimiter=',')
               
            width = 0.35 
            fig, ax = plt.subplots()
            rects1 = ax.bar(np.arange(3), [np.asscalar(data.implicit.values), np.asscalar(data.analytic.values), np.asscalar(data.analytic.values)/np.asscalar(data.implicit.values)], width)#, color='r')
            fontSizes = 20
            ax.set_title('Timing Numerical and analytical', fontsize = fontSizes )
            ax.set_xticks(np.arange(3) + width / 10)
            ax.set_xticklabels(('Implicit', 'Analytical', r'$\frac{Analytical}{Numerical}$'))
            #plt.show()
            fig.tight_layout()
            plt.savefig( outfileName2 + 'Timing.png')
            return data
        
    def project5g(self):
        dx = 0.1#, 0.01]
        safetyFactor = 1.04#[0.96, 1.04]
        dimension = "2D"
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
    
    def project5h(self):
         dxValues = [0.1, 0.02, 0.01, 0.005]#, 0.025, 0.0125, 0.00625, 0.003125]#, 0.01]
         safetyFactor = 0.9#, 1.04]
         movieCounter = 0
         scenario = "2DJacobiIterations"
         threadNumber = 1
         FontSizeUniversal = 22
         theta = 0.5
         T = .3
         data = OrderedDict()
         data2 = OrderedDict()
         counter = 1
         iterationArrayJacobi = []
         iterationArrayGaussSeidel = []
         NValues = []
         for dx in dxValues:
             nx = int(round(1./dx + 1))
             NValues.append(nx)
             dt = dx**2/4.0*(1/safetyFactor)
             alpha = dt/dx**2
             outfileName ='out5h'
             outfileName2 = os.getcwd() + '/results/' + outfileName
             self.runCpp(outfileName2, dt, dx, theta, T, scenario, threadNumber)
             data[counter] = pd.read_csv(outfileName2 + 'JacobiIterationNumber.txt', header=None)# , delimiter=',')
             data2[counter] = pd.read_csv(outfileName2 + 'GaussSeidelIterationNumber.txt', header=None)# , delimiter=',')
             iterationArrayJacobi.append(np.asscalar(data[counter].values))
             iterationArrayGaussSeidel.append(np.asscalar(data2[counter].values))
             counter +=1
         fig, ax = plt.subplots()
         ax.plot(NValues, iterationArrayJacobi, '-o', NValues, iterationArrayGaussSeidel, '-x')
         ax.set_xlabel('N')
         ax.set_ylabel('k/N')
         ax.set_title('Comparison iterative methods')
         ax.legend(['Jacobi', 'Gauss Seidel'], loc=0) 
         fig.tight_layout()
         fig.savefig(outfileName2+ '.png')
         plt.close()
         
         fig, ax = plt.subplots()
         ax.plot(NValues, np.asarray(iterationArrayGaussSeidel)/np.asarray(iterationArrayJacobi), '-o')
         ax.set_xlabel('N')
         ax.set_ylabel(r'$k_{GS}/k_{Jacobi}$')
         ax.set_title('Comparison iterative methods')
         #ax.legend(['Jacobi', 'Gauss Seidel'], loc=0) 
         fig.tight_layout()
         fig.savefig(outfileName2+ '2.png')
         plt.close()
         
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

    if args.task  == '5h':
       project5h = Project5()
       data = project5h.project5h()
       
    if args.task  == '2dnorm':
       project5d2D = Project5()
       data = project5d2D.project5d2D()
       
#%%



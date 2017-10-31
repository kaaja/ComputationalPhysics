#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict
import matplotlib.mlab as mlab


#%% clean results and movie directory
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

#%% Run 
class Project4:
    
    def __init__(self):
        return None
        
    def runCpp(self, outfileName, n_spins,  mcs, initial_temp,
         final_temp, temp_step, orderingFixed):
        """
        Compiles and runs cpp program from the command line and
        makes sure the mesh size is not too big.
        """
        n_spins = str(n_spins)
        mcs = str(mcs)
        initial_temp = str(initial_temp )
        final_temp = str(final_temp) 
        temp_step = str(temp_step)
        call(["./AllrunVectorized", outfileName, n_spins,  mcs, initial_temp,
         final_temp, temp_step, orderingFixed])
    
    def project4b(self):
        #results4b =  OrderedDict()
        N = 5 # Number of different sizes for the MC experiments
        MCSamples = [10**i for i in xrange(2,N)]
        outfileName = 'Out4b'
        n_spins = 2
        initial_temp = 1.0
        final_temp = 1.
        temp_step = 0.1
        orderingFixed = 'orderingFixed'
        
        for mcs in MCSamples:
            outfileName2 = os.getcwd() + '/results/' + outfileName
            outfileName3 = outfileName2 + 'Temp' + str(initial_temp).replace(".", "").replace("0", "") + '.csv'
            self.runCpp(outfileName2, n_spins,  mcs, initial_temp,
             final_temp, temp_step, orderingFixed)
        results4b = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
        results4b.columns = ["acceptedMoves","mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
        #results4b.to_csv(outfileName2 + '4Rport', sep='\t', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        results4b.to_latex(outfileName2 + '4Rport.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)

        results4b2 = results4b.copy()
        
        Z = 4.*(3.+np.cosh(8.))
        EavgExact = -32.*np.sinh(8.)/Z/n_spins/n_spins
        absMavgExact = 8.*(np.exp(8.)+2.)/Z/n_spins/n_spins
        CvExact = (256.*np.cosh(8.)-32.**2*(np.sinh(8.))**2./Z)/Z/n_spins/n_spins
        chiExact = 32.*(np.exp(8.)+1.-2.*(np.exp(8.)+2.)**2./Z)/(Z*initial_temp)/n_spins/n_spins
        
        results4b2.Eavg = (results4b2.Eavg/EavgExact-1.)*100
        results4b2.absMavg = (results4b2.absMavg/absMavgExact-1.)*100
        results4b2.Cv = (results4b2.Cv/CvExact-1.)*100
        results4b2.chi = (results4b2.chi/chiExact-1.)*100
        
        #results4b2.to_csv(outfileName2 + '4RportDev', sep='\t', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        results4b2.to_latex(outfileName2 + '4RportDev.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        return results4b, results4b2
    
    def project4c(self):
        results4cFixed=  OrderedDict()
        results4cRandom  = OrderedDict()
        N = 10 # Number of different sizes for the MC experiments
        MCSamples = [100*2**i for i in xrange(0,N)]
        outfileName = 'Out4c'
        n_spins = 20
        initial_temp = 1.
        final_temp = 2.4
        temp_step = 1.4
        numberOfTemperatures = 2
        orderingTypes = ['orderingFixed', 'nonfixed']
        temperatures = [initial_temp + i*temp_step for i in xrange(numberOfTemperatures)]
        
        
        for orderingType in orderingTypes:
            for mcs in MCSamples:
                outfileName2 = os.getcwd() + '/results/' + orderingType + outfileName 
                self.runCpp(outfileName2, n_spins,  mcs, initial_temp,
                 final_temp, temp_step, orderingType)
            if orderingType == 'orderingFixed':
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + '.csv'    
                    results4cFixed['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4cFixed['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    results4cFixed['temperature %f' % temperature].to_latex(outfileName2 + '4Rport.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
            else:
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + '.csv'    
                    results4cRandom['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4cRandom['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    results4cRandom['temperature %f' % temperature].to_latex(outfileName2 + '4Rport.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)

        fig4, ax4 = plt.subplots()
        ax4.hold('on')
        ax4.set_xlabel(r'$\log_{2} MCS$')
        ax4.set_ylabel(r"$\log_{2}\frac{Accepted moves}{L^2}$")
        legends4 = []
        
        for temperature in temperatures:
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            
            fig5, ax5 = plt.subplots()

            ax2.set_xlabel(r'$\log_{2} MCS$')
            ax2.set_ylabel(r"$\frac{<E>}{l^2}$")
            ax3.set_xlabel(r'$\log_{2} MCS$')
            ax3.set_ylabel(r"$\frac{<|M|>}{l^2}$")
            
            ax5.set_xlabel(r'$\log_{2} MCS$')
            ax5.set_ylabel(r"$\log_{2}\frac{Accepted moves}{L^2}$")

            legends = ['Fixed', 'Random']
            
            acceptedMovesShareFixed = results4cFixed['temperature %f' % temperature].acceptedMoves/float(n_spins)**2 
            acceptedMovesShareRandom = results4cRandom['temperature %f' % temperature].acceptedMoves/float(n_spins)**2
            
            ax2.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].Eavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].Eavg)
            ax3.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].absMavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].absMavg)
            
            ax5.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), np.log2(acceptedMovesShareFixed) , 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), np.log2(acceptedMovesShareRandom))
            
            
            ax4.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), np.log2(acceptedMovesShareFixed) , 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), np.log2(acceptedMovesShareRandom))
            legends4.append('Initial fixed,  temp = %.2f' %temperature)
            legends4.append('Initial random, temp = %.2f' %temperature)
            
            ax2.legend(legends, loc=0)
            ax3.legend(legends, loc=0)
            ax5.legend(legends, loc=0)
            
            ax2.set_title('Temperature = %.2f' %temperature)
            ax3.set_title('Temperature = %.2f' %temperature)
            ax5.set_title('Temperature = %.2f' %temperature)
            ax2.get_yaxis().get_major_formatter().set_useOffset(False)
            ax3.get_yaxis().get_major_formatter().set_useOffset(False)
            ax5.get_yaxis().get_major_formatter().set_useOffset(False)
            fig2.tight_layout()
            fig3.tight_layout()
            fig5.tight_layout()
            fig2.savefig('results/4cEnergy' + str(temperature).replace(".", "").replace("0", "") + '.png') 
            fig3.savefig('results/4cMoment'  + str(temperature).replace(".", "").replace("0", "") + '.png') 
            fig5.savefig('results/4cAcceptedMoves'  + str(temperature).replace(".", "").replace("0", "") + '.png') 
        
        ax4.legend(legends4, loc=0)
        ax4.get_yaxis().get_major_formatter().set_useOffset(False)
        fig4.tight_layout()
        fig4.savefig('results/4cAcceptedMoves.png') 
        plt.close()
        

        return results4cFixed, results4cRandom

    def project4d(self):
        results4dFixed=  OrderedDict()
        results4dRandom  = OrderedDict()
        results4dFixedEnergyArray = OrderedDict()
        results4dRandomEnergyArray = OrderedDict()
        N = 16 # Number of different sizes for the MC experiments
        MCSamples = [2**N]
        outfileName = 'Out4d'
        n_spins = 20
        initial_temp = 1.
        final_temp = 2.4
        temp_step = 1.4
        numberOfTemperatures = 2
        orderingTypes = ['orderingFixed', 'nonfixed']
        temperatures = [initial_temp + i*temp_step for i in xrange(numberOfTemperatures)]
        
        
        for orderingType in orderingTypes:
            for mcs in MCSamples:
                outfileName2 = os.getcwd() + '/results/' + orderingType + outfileName 
                self.runCpp(outfileName2, n_spins,  mcs, initial_temp,
                 final_temp, temp_step, orderingType)
            if orderingType == 'orderingFixed':
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") +"Mcs%d" %mcs + '.csv'    
                    results4dFixedEnergyArray['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4dFixedEnergyArray['temperature %f' % temperature].columns = ["energy"]
                    
                    outfileName4 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + '.csv'    
                    results4dFixed['temperature %f' % temperature] = pd.read_csv(outfileName4, delim_whitespace=True, header=None)
                    results4dFixed['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    
            else:
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + "Mcs%d" %mcs +'.csv'    
                    results4dRandomEnergyArray['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4dRandomEnergyArray['temperature %f' % temperature].columns = ["energy"]
                    
                    outfileName4 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + '.csv'    
                    results4dRandom['temperature %f' % temperature] = pd.read_csv(outfileName4, delim_whitespace=True, header=None)
                    results4dRandom['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]


        

        
        
        for temperature in temperatures:
            fig6, ax6 = plt.subplots()
            #ax6.hold('on')
            ax6.set_xlabel(r'$Energy$')
            ax6.set_ylabel(r"$Probability$")
            
            fig2, ax2 = plt.subplots()
            #ax6.hold('on')
            ax2.set_xlabel(r'$Energy$')
            ax2.set_ylabel(r"$Probability$")

#            ax6.hist(results4dFixed['temperature %f' % temperature].energy.values)#, density=True, stacked=True)
#            ax2.hist(results4dRandom['temperature %f' % temperature].energy.values)#, density=True, stacked=True)
                    
            weights = np.ones_like(results4dFixedEnergyArray['temperature %f' % temperature].energy.values)/float(len(results4dFixedEnergyArray['temperature %f' % temperature].energy.values))
            ax6.hist(results4dFixedEnergyArray['temperature %f' % temperature].energy.values, bins = 25, weights=weights, edgecolor='k')
            weights = np.ones_like(results4dRandomEnergyArray['temperature %f' % temperature].energy.values)/float(len(results4dRandomEnergyArray['temperature %f' % temperature].energy.values))
            #ax2.hist(results4dRandom['temperature %f' % temperature].energy.values, bins = 25, weights=weights, edgecolor='k')
            #ax2 = seaborn.distplot(results4dRandom['temperature %f' % temperature].energy.values, hist_kws=dict(edgecolor="k", linewidth=2), norm_hist=True)
            num_bins = 25
            n, bins, patches = ax2.hist(results4dRandomEnergyArray['temperature %f' % temperature].energy.values, num_bins, normed=1,  edgecolor='k')#, weights = weights)
            mu = results4dRandom['temperature 2.400000'].Eavg.values*n_spins**2
            sigma = np.sqrt(results4dRandom['temperature 2.400000'].sigmaE.values*n_spins**2)
            y = mlab.normpdf(bins, mu, sigma)
            ax2.plot(bins, y, '--')

            ax6.set_title('Fixed initial config \n Temperature = %.2f' %temperature)
            ax6.get_yaxis().get_major_formatter().set_useOffset(False)
            fig6.tight_layout()
            
            ax2.set_title(r'Random initial config' '\n' r'Temperature = %.2f' ' \n' r'$\mu = %.2f' ' ' r'  \sigma = %.2f$' %(temperature, mu, sigma))
            ax2.get_yaxis().get_major_formatter().set_useOffset(False)
            fig2.tight_layout()

            fig6.savefig('results/4dHistogramFixed'  + str(temperature).replace(".", "").replace("0", "") + '.png') 
            fig2.savefig('results/4dHistogramRandom'  + str(temperature).replace(".", "").replace("0", "") + '.png') 

            #plt.close()
        plt.close()
        

        return results4dFixed, results4dRandom


        
        
#%%
scenario = '4d'

if scenario == '4b':
    project4b = Project4()
    results4b, results4b2 = project4b.project4b()
elif scenario == '4c':
    project4c = Project4()
    results4cFixed, results4cRandom = project4c.project4c()
elif scenario == '4d':
    project4d = Project4()
    results4dFixed, results4dRandom = project4d.project4d()





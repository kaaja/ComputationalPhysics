#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict


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
        N = 8 # Number of different sizes for the MC experiments
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
        N = 5 # Number of different sizes for the MC experiments
        MCSamples = [10**i for i in xrange(2,N)]
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
        ax4.set_ylabel(r"$\frac{Accepted moves}{mcs /cdotl^2} \cdot 100.$")
        legends4 = []
        
        for temperature in temperatures:
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()

            ax2.set_xlabel(r'$\log_{2} MCS$')
            ax2.set_ylabel(r"$\frac{<E>}{l^2}$")
            ax3.set_xlabel(r'$\log_{2} MCS$')
            ax3.set_ylabel(r"$\frac{<|M|>}{l^2}$")
            
            legends = ['Fixed', 'Random']
            
            ax2.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].Eavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].Eavg)
            ax3.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].absMavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].absMavg)
            
            
            acceptedMovesShareFixed = results4cFixed['temperature %f' % temperature].acceptedMoves/(results4cFixed['temperature %f' % temperature].mcs * float(n_spins)**2 )*100
            acceptedMovesShareRandom = results4cRandom['temperature %f' % temperature].acceptedMoves/(results4cRandom['temperature %f' % temperature].mcs * float(n_spins)**2 )*100
            ax4.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), acceptedMovesShareFixed , 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), acceptedMovesShareRandom )
            legends4.append('Initial fixed,  temp = %.2f' %temperature)
            legends4.append('Initial random, temp = %.2f' %temperature)
            
            ax2.legend(legends, loc=0)
            ax3.legend(legends, loc=0)
            ax2.set_title('Temperature = %.2f' %temperature)
            ax3.set_title('Temperature = %.2f' %temperature)
            ax2.get_yaxis().get_major_formatter().set_useOffset(False)
            ax3.get_yaxis().get_major_formatter().set_useOffset(False)
            fig2.tight_layout()
            fig3.tight_layout()
            fig2.savefig('results/4cEnergy' + str(temperature).replace(".", "").replace("0", "") + '.png') 
            fig3.savefig('results/4cMoment'  + str(temperature).replace(".", "").replace("0", "") + '.png') 
        
        ax4.legend(legends4, loc=0)
        ax4.get_yaxis().get_major_formatter().set_useOffset(False)
        fig4.tight_layout()
        fig4.savefig('results/4cAcceptedMoves.png') 
        plt.close()
        

        return results4cFixed, results4cRandom


        
        
#%%
scenario = '4c'

if scenario == '4b':
    project4b = Project4()
    results4b, results4b2 = project4b.project4b()
elif scenario == '4c':
    project4c = Project4()
    results4cFixed, results4cRandom = project4c.project4c()





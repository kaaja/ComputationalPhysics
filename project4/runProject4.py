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
         final_temp, temp_step):
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
         final_temp, temp_step])
    
    def project4b(self):
        #results4b =  OrderedDict()
        N = 11 # Number of different sizes for the MC experiments
        MCSamples = [10**i for i in xrange(2,N)]
        outfileName = 'Out4b.csv'
        n_spins = 2
        initial_temp = 1.
        final_temp = 1.
        temp_step = 0.1
        
        for mcs in MCSamples:
            outfileName2 = os.getcwd() + '/results/' + outfileName
            self.runCpp(outfileName2, n_spins,  mcs, initial_temp,
             final_temp, temp_step)
        results4b = pd.read_csv(outfileName2, delim_whitespace=True, header=None)
        results4b.columns = ["mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
        #results4b.to_csv(outfileName2 + '4Rport', sep='\t', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        results4b.to_latex(outfileName2 + '4Rport', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)

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
        results4b2.to_latex(outfileName2 + '4RportDev', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        return results4b, results4b2

        
        
#%%
project4b = Project4()
results4b, results4b2 = project4b.project4b()



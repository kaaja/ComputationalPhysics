#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from subprocess import call
import os
from collections import OrderedDict


#%% clean results and movie directory
#call(["./Allclean"]) 
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
        return 0
        
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
        
#%%
project4b = Project4()
outfileName = 'OutPythonTest'
n_spins = 4
mcs = 10000
initial_temp = 1.
final_temp = 1.1
temp_step = 0.1


project4b.runCpp(outfileName=outfileName, n_spins=n_spins,  mcs=mcs, initial_temp=initial_temp,
         final_temp=final_temp, temp_step=temp_step)



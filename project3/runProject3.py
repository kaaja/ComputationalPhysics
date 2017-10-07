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

#%% 2, run 

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
    outfileName = 'sunEarth.csv'
    N = 1000
    finalTime = 1.0
    solverType = 'forwardEuler'
    runCpp(outfileName, finalTime, N, solverType)
    #sunEarth= OrderedDict()
    sunEarth = pd.read_table("results/" + outfileName, 
        			            delimiter=',')
    
    plotSunEarth(sunEarth)
    return sunEarth

def plotSunEarth(sunEarth):
    fig, ax = plt.subplots()
    ax.plot(sunEarth.x, sunEarth.y)
    ax.set_title("Sun Earth. Forward Euler")
    ax.set_xlabel('x [Au]')
    ax.set_ylabel('y [Au]')
    ax.set_xlim(-2., 2.)
    ax.set_ylim(-2., 2.)
    #outfileName = '%sOmega%sRhoMaxComparison' %(electronType,  omega.replace(".", ""))
    filename = ('sunEarth.png')
    fig.savefig(filename) 
    plt.show()
    #plt.close()


#%% 2, run     
sunEarth()

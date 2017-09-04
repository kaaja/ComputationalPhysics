# Importing output-files from "warm_up2.cpp"
# Create figure
# Estimate slopes

#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

#%%

#data = pd.read_csv("/home/karl/doc/subj/att/fys4150/week1_warm_up/out_warm_up2.csv") # Seems to need full addres

gaussianTridiagonalScalars = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_scalars", 
                     delimiter=',') # Seems to need full addres
gaussianTridiagonalSymmetricScalars = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonalSymmetric_scalars", 
                     delimiter=',') # Seems to need full addres
LUScalars = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/luLib_scalars", 
                     delimiter=',') # Seems to need full addres                                                    
                          
gaussianTridiagonalScalars.columns = gaussianTridiagonalScalars.columns.str.strip().str.replace(' ', '_') # Fixing white space issue variable names
gaussianTridiagonalSymmetricScalars.columns = gaussianTridiagonalSymmetricScalars.columns.str.strip().str.replace(' ', '_')
LUScalars.columns = LUScalars.columns.str.strip().str.replace(' ', '_')

gaussianTridiagonalLogTimes = np.log10(gaussianTridiagonalScalars.time_used)
gaussianTridiagonalSymmetricLogTimes = np.log10(gaussianTridiagonalSymmetricScalars.time_used)
LuLogTimes = np.log10(LUScalars.time_used)
#logGaussianTridiangoanlTime = np.log10(gaussianTridiagonalScalars.time_used)

#%% Plot of log times
plt.figure()
plt.plot(gaussianTridiagonalScalars.log_h, gaussianTridiagonalLogTimes)
plt.hold('on')
plt.plot(gaussianTridiagonalSymmetricScalars.log_h, gaussianTridiagonalSymmetricLogTimes)
plt.plot(LUScalars.log_h, LuLogTimes)         
plt.legend(['Thomas', 'Symmetric', 'LU'], fontsize = 'large', loc = 'upper right')
plt.title('CPU times', fontsize = 'xx-large')
plt.xlabel('log h', fontsize = 'xx-large')
plt.ylabel('log time', fontsize = 'xx-large')
plt.savefig('/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/logTimes.pdf')



#%%
numberOfSimulations = len(gaussianTridiagonalScalars)
x = {}
gaussianTridiagonal = {}
gaussianTridiagonalSymmetric = {}
exactSolution = {}

for key in xrange(numberOfSimulations):    
    gaussianTridiagonal[key] = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_%d" %(key+1), 
                     delimiter=',').values
    gaussianTridiagonal[key] = np.reshape(gaussianTridiagonal[key], -1)
    gaussianTridiagonalSymmetric[key] = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonalSymmetric_%d" %(key+1), 
                     delimiter=',').values
    gaussianTridiagonalSymmetric[key] = np.reshape(gaussianTridiagonalSymmetric[key], -1)
    N = len(gaussianTridiagonal[key])
    h = 1./(N+1)
    x[key] = np.linspace(h, 1.-h, N)


                        
#%% Plot cpp-output
plt.figure()
x = np.linspace(0,1,len(exact_2))
plt.plot(x, exact_2, x, numerical_2)
plt.figure()
x2 = np.linspace(0,1,len(exact_1))
plt.plot(x2, exact_1, x2, numerical_1)
                        
#%% Plot cpp-output
plt.plot(data.log_h, data.log_rel_error)
plt.xlabel('log h')
plt.ylabel('log rel error')
plt.show()

#%% Estimation of slopes
data2 = data.values # get rid of pandas indices etx. 

index_minumim_error_h = np.argmin(data2[:,1]) 
xdata_left = data2[0:index_minumim_error_h,0]
ydata_left = data2[0:index_minumim_error_h,1]

f = lambda x, a, b: a*x + b # Assume linear model

a0 = .4 # Initiel guess for estimator
b0 = 2.

xdata = xdata_left
ydata = ydata_left
popt, pcov = curve_fit(f, xdata, ydata, p0=[a0, b0])
residuals = ydata - f(xdata,*popt)

plt.figure()
plt.plot(xdata, f(xdata, popt[0], popt[1]), label=('a=%.1f, b= %.4f' %(popt[0], popt[1])))
plt.hold('on')
plt.scatter(xdata, ydata, marker='*', c='b', label='data') #s=1, 
plt.show()

#%% Test of error calculations
absolute_error = numerical_2.values-exact_2.values
plt.plot(x, abs(absolute_error/exact_2.values))
max_error = max(absolute_error)
logError = (np.log10(np.abs(absolute_error/exact_2.values)))
maxLogError = max(np.log10(np.abs(absolute_error/exact_2.values)))

l2 = np.sqrt(np.sum(absolute_error**2)/len(absolute_error))
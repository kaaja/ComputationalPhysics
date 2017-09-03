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
gaussianTridiagonalScalarData = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_scalars", 
                     delimiter=',') # Seems to need full addres
                                            
gaussianSymmetricTridiagonalScalarData = pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonalSymmetric_scalars", 
                     delimiter=',') # Seems to need full addres
#%% Reading of arrays, numerical and exact solutions
gaussianTridiagonal = {}
gaussianTridiagonalSymmetric = {}
exact_solutions = {}
number_of_simulations = len(gaussianTridiagonalScalarData)
x = {}
for key in xrange(number_of_simulations):
    gaussianTridiagonal[key] = (pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_numerical%d" %(key+1), 
                      delim_whitespace=True)).values
    gaussianTridiagonal[key] =  np.reshape(gaussianTridiagonal[key], -1)
    exact_solutions[key] = (pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_exact%d" %(key+1), 
                      delim_whitespace=True)).values
    gaussianTridiagonalSymmetric[key] = (pd.read_table("/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonalSymmetric_numerical%d" %(key+1), 
                      delim_whitespace=True)).values
    gaussianTridiagonalSymmetric[key] =  np.reshape(gaussianTridiagonalSymmetric[key], -1)
                   
    exact_solutions[key] =  np.reshape(exact_solutions[key], -1)
    N = len(gaussianTridiagonal[key])
    h = 1./(N+1)
    x[key] = np.linspace(h, 1.-h, N)


#%% Plotting numerical and exact solution
plt.figure()
labels = []
for i in xrange(3):
    plt.plot(x[i], gaussianTridiagonal[i])
    labels.append('N = %d' %len(gaussianTridiagonal[i]))
plt.plot(x[i], exact_solutions[i])
labels.append('Exact')
#plt.legend(labels)
plt.legend(labels, fontsize = 'x-large', loc = 'upper right')
plt.title(' Tridagonal Gaussian', fontsize = 'x-large')
plt.xlabel('x',fontsize = 'x-large')
plt.ylabel('v(x)',fontsize = 'x-large')
plt.savefig('/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_comparison.pdf')
plt.show()

#%% Plot errors
plt.figure()
plt.plot(gaussianTridiagonalScalarData.log_h, gaussianTridiagonalScalarData.log_rel_error)
plt.title(' Thomas algorithm', fontsize = 'x-large')
plt.xlabel('log h',fontsize = 'x-large')
plt.ylabel('log relative error',fontsize = 'x-large')
plt.ylim((-9,-1))
plt.savefig('/home/karl/doc/subj/att/fys4150/build-project1qt-Desktop_Qt_5_9_1_GCC_64bit-Debug/gaussianTridiagonal_error2.pdf')
plt.show()
                        
                        
#%% Errors

#exact_solution = lambda x: 1. - (1. - np.exp(-10.))*x - np.exp(-10.*x)

errors = {}
x = {}
for key in xrange(10):
#%% Plot cpp-output
    N = len(numerical_solutions[key])
    h = 1./(N+1)
    x[key] = np.linspace(h, 1.-h, N)
    numerical_solutions[key] =  np.reshape(numerical_solutions[key], -1)
    errors[key] = (numerical_solutions[key] - exact_solution(x[key]))
    

plt.plot(x[0], numerical_solutions[0], x[0], exact_solution(x[0]))
print numerical_solutions[0], exact_solution(x[0])

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

#%% Comparison different algorithms
x = np.linspace(0,1,len(exact_2))
plt.figure()
plt.subplot(2,1,1)
plt.plot(x, exact_1, x, exact_2)
plt.subplot(2,1,2)
plt.plot(x, (exact_1-exact_2)/exact_1)
(exact_1-exact_2)/exact_1
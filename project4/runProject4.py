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
class Project4:
    
    def __init__(self):
        return None
        
    def runCpp(self,  numprocs, outfileName, n_spins,  mcs, initial_temp,
         final_temp, temp_step, orderingFixed, initializeForEachTemperature, 
         printEnergyArray):
        """
        Compiles and runs cpp program from the command line and
        makes sure the mesh size is not too big.
        """
        n_spins = str(n_spins)
        mcs = str(mcs)
        initial_temp = str(initial_temp )
        final_temp = str(final_temp) 
        temp_step = str(temp_step)
        numprocs = str(numprocs)
        initializeForEachTemperature = str(initializeForEachTemperature)
        call(["./AllrunVectorized", numprocs, outfileName, n_spins,  mcs, initial_temp,
         final_temp, temp_step, orderingFixed, initializeForEachTemperature, 
         printEnergyArray])
    
    def project4b(self):
        N =9 # Number of different sizes for the MC experiments
        MCSamples = [10**i for i in xrange(2,N)]
        outfileName = 'Out4b'
        n_spins = 2
        initial_temp = 1.0
        final_temp = 1.
        temp_step = 0.1
        orderingFixed = 'orderingFixed'
        numprocs = 4
        initializeForEachTemperature = 'notInitializeForEachTemperature'
        printEnergyArray = 'NotprintEnergyArray'
        
        counter = 0
        
        for mcs in MCSamples:
            outfileName2 = os.getcwd() + '/results/' + outfileName
            outfileName3 = outfileName2 + 'TempNumber' + str(counter)+'.csv'
            self.runCpp(numprocs,outfileName2, n_spins,  mcs, initial_temp,
             final_temp, temp_step, orderingFixed, initializeForEachTemperature,
             printEnergyArray)
        results4b = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
        results4b.columns = ["acceptedMoves","mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
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
        
        results4b2.to_latex(outfileName2 + '4RportDev.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
        
        return results4b, results4b2
    
    def project4c(self):
        results4cFixed=  OrderedDict()
        results4cRandom  = OrderedDict()
        N = 22 # Number of different sizes for the MC experiments
        MCSamples = [2**i for i in xrange(7,N)]
        outfileName = 'Out4c'
        n_spins = 20
        initial_temp = 1.
        final_temp = 2.4
        temp_step = 1.4
        numberOfTemperatures = 2
        orderingTypes = ['orderingFixed', 'nonfixed']
        temperatures = [initial_temp + i*temp_step for i in xrange(numberOfTemperatures)]
        numprocs = 4
        initializeForEachTemperature = 'initializeForEachTemperature'
        printEnergyArray = 'NotprintEnergyArray'
        
        for orderingType in orderingTypes:
            for mcs in MCSamples:
                outfileName2 = os.getcwd() + '/results/' + orderingType + outfileName 
                self.runCpp(numprocs,outfileName2, n_spins,  mcs, initial_temp,
                 final_temp, temp_step, orderingType, initializeForEachTemperature, printEnergyArray)
            counter = 0
            if orderingType == 'orderingFixed':
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'TempNumber' + str(counter) + '.csv'    
                    results4cFixed['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4cFixed['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    results4cFixed['temperature %f' % temperature].to_latex(outfileName2 + '4Rport.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
                    counter += 1
            else:
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'TempNumber' + str(counter) + '.csv'    
                    results4cRandom['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4cRandom['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    results4cRandom['temperature %f' % temperature].to_latex(outfileName2 + '4Rport.txt', columns = ['mcs','Eavg', 'absMavg', 'Cv', 'chi'], index=False)
                    counter += 1
                    
        fig4, ax4 = plt.subplots()
        ax4.hold('on')
        ax4.set_xlabel(r'$\log_{2} MCS$')
        ax4.set_ylabel(r"$\frac{Accepted moves}{MC cycles \cdot L^2} \cdot 100$")
        legends4 = []
        
        for temperature in temperatures:
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            
            fig5, ax5 = plt.subplots()

            ax2.set_xlabel(r'$\log_{2} MCS$')
            ax2.set_ylabel(r"$\frac{<E>}{L^2}$")
            ax3.set_xlabel(r'$\log_{2} MCS$')
            ax3.set_ylabel(r"$\frac{<|M|>}{L^2}$")
            
            ax5.set_xlabel(r'$\log_{2} MCS$')
            ax5.set_ylabel(r"$\frac{Accepted moves}{MC cycles \cdot L^2} \cdot 100$")

            legends = ['Fixed', 'Random']
            
            data = results4cFixed['temperature %f' % temperature]
            acceptedMovesShareFixed = [data.acceptedMoves[i]/(float(n_spins)**2*data.mcs[i])*100 for i in  xrange(len(MCSamples))]
            data = results4cRandom['temperature %f' % temperature]
            acceptedMovesShareRandom = [data.acceptedMoves[i]/(float(n_spins)**2*data.mcs[i])*100 for i in  xrange(len(MCSamples))]
            
            ax2.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].Eavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].Eavg)
            ax3.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), results4cFixed['temperature %f' % temperature].absMavg, 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), results4cRandom['temperature %f' % temperature].absMavg)
            
            ax5.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), acceptedMovesShareFixed , 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), acceptedMovesShareRandom)
            
            
            ax4.plot(np.log2(results4cFixed['temperature %f' % temperature].mcs), acceptedMovesShareFixed , 
                    np.log2(results4cRandom['temperature %f' % temperature].mcs), acceptedMovesShareRandom)
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

    def project4d(self, numprocs):
        results4dFixed=  OrderedDict()
        results4dRandom  = OrderedDict()
        results4dFixedEnergyArray = OrderedDict()
        results4dRandomEnergyArray = OrderedDict()
        N = 20 # Number of different sizes for the MC experiments
        MCSamples = [2**N]
        outfileName = 'Out4d'
        n_spins = 20
        initial_temp = 1.
        final_temp = 2.4
        temp_step = 1.4
        numberOfTemperatures = 2
        orderingTypes = ['orderingFixed', 'nonfixed']
        temperatures = [initial_temp + i*temp_step for i in xrange(numberOfTemperatures)]
        numprocs = 4
        initializeForEachTemperature = 'initializeForEachTemperature'
        printEnergyArray = 'printEnergyArray'
        
        for orderingType in orderingTypes:
            for mcs in MCSamples:
                outfileName2 = os.getcwd() + '/results/' + orderingType + outfileName 
                self.runCpp(numprocs, outfileName2, n_spins,  mcs, initial_temp,
                 final_temp, temp_step, orderingType, initializeForEachTemperature,
                 printEnergyArray)
            counter = 0
            if orderingType == 'orderingFixed':
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") +"Mcs%d" %mcs + '.csv'    
                    results4dFixedEnergyArray['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4dFixedEnergyArray['temperature %f' % temperature].columns = ["energy"]
                    
                    outfileName4 = outfileName2 + 'TempNumber' + str(counter)+'.csv'    
                    results4dFixed['temperature %f' % temperature] = pd.read_csv(outfileName4, delim_whitespace=True, header=None)
                    results4dFixed['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    counter += 1

            else:
                for temperature in temperatures:
                    outfileName3 = outfileName2 + 'Temp' + str(temperature).replace(".", "").replace("0", "") + "Mcs%d" %mcs +'.csv'    
                    results4dRandomEnergyArray['temperature %f' % temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                    results4dRandomEnergyArray['temperature %f' % temperature].columns = ["energy"]
                    
                    outfileName4 = outfileName2 + 'TempNumber' + str(counter)+ '.csv'    
                    results4dRandom['temperature %f' % temperature] = pd.read_csv(outfileName4, delim_whitespace=True, header=None)
                    results4dRandom['temperature %f' % temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                    counter += 1

        tableDict = OrderedDict()
        tableDict['Fixed'] = OrderedDict()
        tableDict['Random'] = OrderedDict()
        for temperature in temperatures:
            fig6, ax6 = plt.subplots()
            ax6.set_xlabel(r'$Energy/L^2$')
            ax6.set_ylabel(r"$Probability$")
            
            fig2, ax2 = plt.subplots()
            ax2.set_xlabel(r'$Energy/L^2$')
            ax2.set_ylabel(r"$Probability$")


            
            weights = np.ones_like(results4dRandomEnergyArray['temperature %f' % temperature].energy.values)/float(len(results4dRandomEnergyArray['temperature %f' % temperature].energy.values))
            #ax2 = seaborn.distplot(results4dRandom['temperature %f' % temperature].energy.values, hist_kws=dict(edgecolor="k", linewidth=2), norm_hist=True)
            ax2.hist(results4dRandomEnergyArray['temperature %f' % temperature].energy.values/n_spins**2, bins = 200, range = (-2, -0.5), weights=weights, edgecolor='k') #bins = 25
            #num_bins = 25
            #n, bins, patches = ax2.hist(results4dRandomEnergyArray['temperature %f' % temperature].energy.values/n_spins**2, num_bins, normed=1,  edgecolor='k')#, weights = weights)
            muCpp = np.asscalar(results4dRandom['temperature %f' %temperature].Eavg.values)
            sigmaCpp = np.asscalar(results4dRandom['temperature %f' %temperature].sigmaE.values)
            absMAvg = np.asscalar(results4dRandom['temperature %f' %temperature].absMavg.values)
            sigmaM = np.asscalar(results4dRandom['temperature %f' %temperature].sigmaM.values)
            Cv = np.asscalar(results4dRandom['temperature %f' %temperature].Cv.values)
            chi = np.asscalar(results4dRandom['temperature %f' %temperature].chi.values)
            data = results4dRandomEnergyArray['temperature %f' % temperature].energy.values            
            mu, sigma = norm.fit(data/n_spins**2)            
#            if temperature != 1.00:
#                y = mlab.normpdf(bins, mu, sigma)
#                ax2.plot(bins, y, '--')
            
            #sigma = sigma**2/n_spins**2
            mu = np.average(data)/n_spins**2
            sigma = np.var(data)/n_spins**2
            muDiff = np.asscalar((mu/muCpp-1.)*100)
            sigDiff = np.asscalar((sigma/sigmaCpp - 1.)*100)
            
            tableDict['Random']['temperature %f' %temperature] = [np.log2(MCSamples[0]), temperature, mu, muCpp, muDiff, sigma, sigmaCpp, sigDiff, absMAvg, sigmaM, Cv, chi ]

            
            weights = np.ones_like(results4dFixedEnergyArray['temperature %f' % temperature].energy.values)/float(len(results4dFixedEnergyArray['temperature %f' % temperature].energy.values))
            ax6.hist(results4dFixedEnergyArray['temperature %f' % temperature].energy.values/n_spins**2, bins = 200, range = (-2, -0.5), weights=weights, edgecolor='k') #bins = 'auto', 
            
            muCpp = np.asscalar(results4dFixed['temperature %f' %temperature].Eavg.values)
            sigmaCpp = np.asscalar(results4dFixed['temperature %f' %temperature].sigmaE.values)
            absMAvg = np.asscalar(results4dFixed['temperature %f' %temperature].absMavg.values)
            sigmaM = np.asscalar(results4dFixed['temperature %f' %temperature].sigmaM.values)
            Cv = np.asscalar(results4dFixed['temperature %f' %temperature].Cv.values)
            chi = np.asscalar(results4dFixed['temperature %f' %temperature].chi.values)
            data = results4dFixedEnergyArray['temperature %f' % temperature].energy.values            
            mu, sigma = norm.fit(data/n_spins**2)            

            mu = np.average(data)/n_spins**2
            sigma = np.var(data)/n_spins**2
            muDiff = np.asscalar((mu/muCpp-1.)*100)
            sigDiff = np.asscalar((sigma/sigmaCpp - 1.)*100)
            
            tableDict['Fixed']['temperature %f' %temperature] = [np.log2(MCSamples[0]), temperature, mu, muCpp, muDiff, sigma, sigmaCpp, sigDiff, absMAvg, sigmaM, Cv, chi ]
            
            ax6.set_title('Fixed initial config \n Temperature = %.2f' %temperature)
            ax6.get_yaxis().get_major_formatter().set_useOffset(False)
            fig6.tight_layout()
            
            ax2.set_title(r'Random initial config' '\n' r'Temperature = %.2f' %temperature)
            ax2.get_yaxis().get_major_formatter().set_useOffset(False)
            fig2.tight_layout()

            fig6.savefig('results/4dHistogramFixed'  + str(temperature).replace(".", "").replace("0", "") + '.png') 
            fig2.savefig('results/4dHistogramRandom'  + str(temperature).replace(".", "").replace("0", "") + '.png') 

        plt.close()
        dfFixed =pd.DataFrame.from_dict(tableDict['Fixed'],orient='index')#, index=tableDict['Fixed'])#.T
        outfileName = os.getcwd() + '/results/' + '4dTableFixed.txt'
        dfFixed.columns = [r"$\log_2 MCs$", "T", r"$\mu_E/L^2$", r"$<E>/L^2$", r"$(\frac{\mu_E/L^2}{<E>/L^2} - 1) \cdot 100$", r"$\sigma_E^2/L^2$", r"$\frac{<E^2> - <E>^2}{L^2}$", r"$(\frac{\sigma_E/L^2}{1/L^2(<E^2> - <E>^2)} - 1)\cdot 100$", r"$<|M|>/L^2$", r"$\frac{<|M|^2>/L^2 - <|M|>^2}{L^2}$", r"$Cv/L^2$", r"$\chi/L^2$"]
        dfFixed.to_latex(outfileName, index=False, escape=False)

        dfRandom =pd.DataFrame.from_dict(tableDict['Random'],orient='index')#, index=tableDict['Random'])#.T
        
        outfileName = os.getcwd() + '/results/' + '4dTableRandom.txt'
        dfRandom.columns = [r"$\log_2 MCs$", "T", r"$\mu_E/L^2$", r"$<E>/L^2$", r"$(\frac{\mu_E/L^2}{<E>/L^2} - 1) \cdot 100$", r"$\sigma_E^2/L^2$", r"$\frac{<E^2> - <E>^2}{L^2}$", r"$(\frac{\sigma_E/L^2}{1/L^2(<E^2> - <E>^2)} - 1)\cdot 100$", r"$<|M|>/L^2$", r"$\frac{<|M|^2>/L^2 - <|M|>^2}{L^2}$", r"$Cv/L^2$", r"$\chi/L^2$"]
        dfRandom.to_latex(outfileName, index=False, escape=False)

        return results4dFixed, results4dRandom, results4dFixedEnergyArray, results4dRandomEnergyArray, tableDict
    
    def project4e(self):
        results4e=  OrderedDict()
        N = 1 # Number of different sizes for the MC experiments
        MCSamples = [int(round(10E6))]#[int(round(10E6))]
        mcs = MCSamples[0]
        outfileName = 'Out4e'
        n_spins_list = [40,60,80,100]#, 80]#,100]#, 140]
        initial_temp = 2.2
        final_temp = 2.3
        temp_step = .0125
        numprocs = 4
        orderingType = 'random'
        numberOfTemperatures = 1 + int(round((final_temp - initial_temp)/temp_step))
        temperatures = [initial_temp + i*temp_step for i in xrange(numberOfTemperatures)]
        initializeForEachTemperature = 'notInitializeForEachTemperature'
        printEnergyArray = 'notPrintEnergyArray'

        for n_spins in n_spins_list: 
            results4e['n_spins %d' % n_spins] = OrderedDict()
            outfileName2 = os.getcwd() + '/results/' + outfileName + str(n_spins)
#            self.runCpp(numprocs, outfileName2, n_spins,  mcs, initial_temp,
#                 final_temp, temp_step, orderingType, initializeForEachTemperature,
#                 printEnergyArray)
            counter = 0
            for temperature in temperatures:
                outfileName3 = outfileName2 + 'TempNumber' + str(counter)+'.csv'    #.replace("0", "")
                results4e['n_spins %d' % n_spins]['temperature %f' %temperature] = pd.read_csv(outfileName3, delim_whitespace=True, header=None)
                results4e['n_spins %d' % n_spins]['temperature %f' %temperature].columns = ["acceptedMoves", "mcs", "temperature", "Eavg", "sigmaE", "Mavg", "sigmaM", "absMavg", "Cv", "chi"]
                counter += 1 # Kja: This was mising...So I guess the dict would not be correct.  
                
        # Plotting
        plotVariables = ['Eavg', 'absMavg', 'Cv', 'chi']
        titles = ["$<E>/L^2$", "$<|M|>/L^2$", r"$C_v/L^2$", r"$\chi/L^2$"]
        variableNumber = 0
        for plotvariable in plotVariables:
            fig1, ax1 = plt.subplots()
            legends1 = []
            for n_spins in n_spins_list:
                data = []
                for temperature in temperatures:
                    if plotvariable == 'Eavg':
                        data.append(results4e['n_spins %d' %n_spins]['temperature %f' %temperature].Eavg)
                    elif plotvariable == 'absMavg':
                        data.append(results4e['n_spins %d' %n_spins]['temperature %f' %temperature].absMavg)
                    elif plotvariable == 'Cv':
                        data.append(results4e['n_spins %d' %n_spins]['temperature %f' %temperature].Cv)
                    elif plotvariable == 'chi':
                        data.append(results4e['n_spins %d' %n_spins]['temperature %f' %temperature].chi)

                ax1.plot(temperatures, data, '-o')
                legends1.append('L %d' %n_spins)
            ax1.set_xlabel(r'$Temperature$')
            ax1.set_ylabel(titles[variableNumber])
            ax1.set_title(titles[variableNumber])
            variableNumber += 1
            ax1.legend(legends1, loc=0)
            fig1.tight_layout()
            fig1.savefig('results/4e'  + plotvariable + '.png') 
            plt.close()
            
        # Critical temperature
        n_spins_list = [40, 60, 80, 100] # new for table only
        TCritical = []
        for n_spins in n_spins_list:
            CvMax = [np.asscalar(results4e['n_spins %d' %n_spins]['temperature %f' %temperature].Cv.values) for temperature in temperatures]
            TCriticalIndex = np.argmax(CvMax)
            TCritical.append(temperatures[TCriticalIndex])
        
        numberOfSpins = len(n_spins_list)
        TCriticalInfList = []
        for numberOfSpinsCriticalTemperature in xrange(2, numberOfSpins+1):
            numberOfCombinations = int(round(comb(numberOfSpinsCriticalTemperature, 2)))
            TCriticalInf = 0
            for i in xrange(numberOfSpinsCriticalTemperature-1):
                j = i+1
                while  j < numberOfSpinsCriticalTemperature:
                    aOverL = (TCritical[i] - TCritical[j])/(1. - n_spins_list[i]/n_spins_list[j])
                    TCriticalInf += TCritical[i] - aOverL
                    j += 1
            
            TCriticalInf = TCriticalInf/numberOfCombinations
            TCriticalInfList.append(TCriticalInf)            
            
        outfileName2 = os.getcwd() + '/results/' + outfileName + 'TcCritical.txt'
        TOnsager = 2./np.log(1 + np.sqrt(2))
        diffFromOnsager = (np.asarray(TCriticalInfList)/TOnsager - 1.)*100
        spinCombos = np.asarray([r'$[40, 60]$', r'$[40, 60, 80]$', r'$[40, 60, 80, 100]$'])
        
        tableForTcCritical = np.transpose(np.array((spinCombos, np.asarray(TCriticalInfList),diffFromOnsager)))
        tableForTcCriticalToLatex  = pd.DataFrame(tableForTcCritical, columns=['Spin combos', r'$T_c^{Estimate}(L=\infty)$',r'$(\frac{T_c^{Estimate}(L=\infty)}{T_{c,exact}}-1) \cdot 100$' ] )
        tableForTcCriticalToLatex.to_latex(outfileName2, index=False, escape=False)
        
        return results4e      
    
    def projectMPITiming(self):
        results4e=  OrderedDict()
        N = 20 # Number of different sizes for the MC experiments
        MCSamples = [2**N]
        
        mcs = MCSamples[0]
        outfileName = 'Time'
        
        n_spins = 20
        numProcsList = [1,2,3,4]
        initial_temp = 1.0
        final_temp = 1.0
        temp_step = .0125
        #numprocs = 4
        orderingType = 'random'
        initializeForEachTemperature = 'notInitializeForEachTemperature'
        printEnergyArray = 'notPrintEnergyArray'

        for numprocs in numProcsList: 
            #results4e['numprocs  %d' % numprocs ] = OrderedDict()
            outfileName2 = os.getcwd() + '/results/' + outfileName 
            self.runCpp(numprocs, outfileName2, n_spins,  mcs, initial_temp,
                 final_temp, temp_step, orderingType, initializeForEachTemperature,
                 printEnergyArray)
            
            outfileName3 = outfileName2 + 'NumProcs' + str(numprocs)+'.csv'    #.replace("0", "")
            results4e['numprocs %d' % numprocs]= pd.read_csv(outfileName3, delim_whitespace=True)

        fig, ax = plt.subplots()
        data = []
        for i in xrange(1,5):
            data.append(np.asscalar(results4e['numprocs %d' %i].totalTime.values))
        ax.plot(np.log2(numProcsList), data, '-o')
        ax.set_xlabel(r'$\log_2$ Number of processors')
        ax.set_ylabel(r"Time [s]")
        fig.tight_layout()
        fig.savefig('results/Timing.png') 
        plt.close()
            
        
        
#        tableForTcCritical = np.transpose(np.array((spinCombos, np.asarray(TCriticalInfList),diffFromOnsager)))
#        tableForTcCriticalToLatex  = pd.DataFrame(tableForTcCritical, columns=['Spin combos', r'$T_c^{Estimate}(L=\infty)$',r'$(\frac{T_c^{Estimate}(L=\infty)}{T_{c,exact}}-1) \cdot 100$' ] )
#        tableForTcCriticalToLatex.to_latex(outfileName2, index=False, escape=False)
        
        return results4e      
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
    
    parser = argparse.ArgumentParser(description="starts a c++ program simulating phase transitions with 2D Ising model, generatsplots. \n Alternatives: 4b, 4c, 4d, 4e")
    parser.add_argument("task", type=str, default='4b', help="press 4b, 4c, 4d, or 4e")
    args = parser.parse_args()
     
    if args.task  == '4b':
       project4b = Project4()
       results4b, results4b2 = project4b.project4b()
    elif args.task == '4c':
       project4c = Project4()
       results4cFixed, results4cRandom = project4c.project4c()
    elif args.task  == '4d':
       numprocs = 4
       project4d = Project4()
       results4dFixed, results4dRandom, results4dFixedEnergyArray, results4dRandomEnergyArray, tableDict= project4d.project4d(numprocs)
    elif args.task == '4e':
       project4e = Project4()
       results4e = project4e.project4e()
    elif args.task == 'mpiTiming':
       project4mpiTiming = Project4()
       results4e = project4mpiTiming.projectMPITiming()
       
#%%



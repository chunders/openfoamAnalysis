#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
       _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Tue Jun 12 10:18:09 2018

@author: chrisunderwood
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class select_Gaus_fit():
    def __init__(self, x, y, plottingOnOff, guess):
        self.x = x
        self.y = y
        self.plotting = plottingOnOff
        self.guess = guess
        if self.plotting:
            self.plotInput()
            
        self.fit_gaus_2_to_8()
        
    def gaus2(self, x, *params):
        #Gaussian function
        A = params[0]
        x0 = params[1]
        c = params[2]
        return A*np.exp(-((x-x0)/(c**2))**2)

    def gaus4(self, x, *params):
        #Gaussian function
        A = params[0]
        x0 = params[1]
        c = params[2]
        return A*np.exp(-((self.x-x0)/(c**2))**4)
    
    def gaus6(self, x, *params):
        #Gaussian function
        A = params[0]
        x0 = params[1]
        c = params[2]
        return A*np.exp(-((self.x-x0)/(c**2))**6)
    
    def gaus8(self, x, *params):
        #Gaussian function
        A = params[0]
        x0 = params[1]
        c = params[2]
        return A*np.exp(-((self.x-x0)/(c**2))**8)
    
    def plotInput(self):
        plt.plot(self.x,self.y)
        
    def fitGaus2(self):
        self.popt_G2, self.pcov_G2 = curve_fit(self.gaus2, self.x, self.y , 
                                               p0=self.guess)
        self.fitG2 = self.gaus2(self.x, *self.popt_G2)
        
        if self.plotting:
            print 'Fitting SG2: '
            print self.guess
            print self.popt_G2
            print
            plt.plot(self.x, self.fitG2)
            
    def fitGaus4(self):
        self.popt_G4, self.pcov_G4 = curve_fit(self.gaus4, self.x, self.y , 
                                               p0=self.guess)
        self.fitG4 = self.gaus4(self.x, *self.popt_G4)        
        if self.plotting:
            print 'Fitting SG4: '
            print self.guess
            print self.popt_G4
            print
            plt.plot(self.x, self.fitG4)
       
    def fitGaus6(self):
        self.popt_G6, self.pcov_G6 = curve_fit(self.gaus6, self.x, self.y , 
                                               p0=self.guess)
        self.fitG6 = self.gaus6(self.x, *self.popt_G6)
        if self.plotting:
            print 'Fitting SG6: '
            print self.guess
            print self.popt_G6
            print
            plt.plot(self.x, self.fitG6)    

    def fitGaus8(self):
        
        self.popt_G8, self.pcov_G8 = curve_fit(self.gaus8, self.x, self.y , 
                                               p0=self.guess)
        self.fitG8 = self.gaus8(self.x, *self.popt_G8)
        if self.plotting:
            print 'Fitting SG8: '
            print self.guess
            print self.popt_G8
            print
            plt.plot(self.x, self.fitG8)   
    
    def fit_gaus_2_to_8(self):
        self.fitOptions = ['SG2', 'SG4', 'SG6', 'SG8']
        self.fitPowers = [2, 4, 6, 8]

        self.fitGaus2()
        self.fitGaus4()
        self.fitGaus6()
        self.fitGaus8()
        self.fitParams = [self.popt_G2, self.popt_G4, self.popt_G6, self.popt_G8]
        self.CalcR_2_value()

    def nearposn(self, array,value):
        posn = (abs(array-value)).argmin()
        return posn

        
    def CalcR_2_value(self):
        self.rrValues = []
        for fit in [self.fitG2, self.fitG4, self.fitG6, self.fitG8]:
            
            residuals = self.y - fit
            ss_res = np.sum(residuals**2)
            ss_tot = np.sum((self.y-np.mean(self.y))**2)
            self.rrValues.append( 1 - (ss_res / ss_tot) )
        maxRR = max(self.rrValues)
        self.bestFit = self.nearposn(self.rrValues, maxRR)
        outStr = ''
        CSI="\x1B[31;40m"
        CEND = '\x1B[0m'
        for i in self.rrValues:
            if i == self.rrValues[self.bestFit]:
                outStr += CSI + str(i) + ' ' + CEND
            else:
                outStr += str(i) + ' '
        if self.plotting: print outStr 
        
    def output(self):
        if self.plotting:
            print self.fitOptions[self.bestFit]
            print self.fitParams[self.bestFit]
            print self.fitPowers[self.bestFit]
        return self.fitPowers[self.bestFit] , self.fitParams[self.bestFit]
             





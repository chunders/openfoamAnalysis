#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
       _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
Created on Wed Jul 18 12:41:27 2018

@author: chrisunderwood

Fitting function to the residual from the big gaussian fits
"""
import numpy as np
import matplotlib.pyplot as plt
import select_best_gaussFit_class as bestGaus


def nearposn(array,value):
	#Find array position of value
    posn = (abs(array-value)).argmin()
    return posn


def gaussian2(x, *params):
    #Gaussian function
    A = params[0]
    x0 = params[1]
    c = params[2]
    return A*np.exp(-((x-x0)/(c**2))**2)

def gaussian4(x, *params):
    #Gaussian function
    A = params[0]
    x0 = params[1]
    c = params[2]
    return A*np.exp(-((x-x0)/(c**2))**4)

def gaussian6(x, *params):
    #Gaussian function
    A = params[0]
    x0 = params[1]
    c = params[2]
    return A*np.exp(-((x-x0)/(c**2))**6)

def gaussian8(x, *params):
    #Gaussian function
    A = params[0]
    x0 = params[1]
    c = params[2]
    return A*np.exp(-((x-x0)/(c**2))**8)

class resFitter():
#==============================================================================
# This class takes in the res from the first fit, and then fits it
# with a series of gaussians to improve on the fit.
# Splitting the fitting up into points where it changes around zero
#==============================================================================
    def __init__(self,xRes, Res, shockFront, Graphics, printing):
        # Initialise the class
        self.shock  = shockFront
        self.x_res = xRes
        self.res = Res
        self.displayToConsole = printing
        # steepness of stepfunctions 
        self.steepness = 1e6
        
        #Work out when it changes sign
        self.nodeIndex = []
        for i in range(len(self.res)-1):
            if self.res[i] * self.res[i+1] < 0:
                self.nodeIndex.append(i)
        
        self.nodesInRes = self.x_res[self.nodeIndex] 
        self.graphics = Graphics
        if self.graphics: self.displayInput()
        
        self.textOutput = ''
        
        
    def displayInput(self):
        plt.plot(self.x_res, self.res)
        plt.plot(self.nodesInRes, np.zeros(len(self.nodesInRes)), 's')
        
        #Display the shock
        plt.vlines(self.shock, min(self.res), max(self.res))
        plt.show()

    def sections(self):
#==============================================================================
#         Need to create a fit for each of the different areas, 
#         the first attempt is to be done with guassians
#==============================================================================
        self.gauss_Fit_params = []
        if self.displayToConsole: print 'Length of the array' , len(self.x_res)
        self.newFit = np.zeros(len(self.x_res))
        self.newX = []
        startIndex = 1e9
        endIndex  = 0
        self.regionCounter = 0
        
#==============================================================================
#         Look at the different regions between nodes
#==============================================================================
        for i in range(len(self.nodeIndex) - 1):
#==============================================================================
#             Only the nodes around the shock are really important
#             Therefore limit the fitting to the region that is important
#==============================================================================
            fitBeyondShock = 5e-05 # In meters
            if self.nodesInRes[i+1] > 0 and self.nodesInRes[i] < self.shock + fitBeyondShock:
                # Count the regions for labelling sake
                self.regionCounter +=1
                if self.nodeIndex[i] < startIndex: startIndex = self.nodeIndex[i]
                if self.nodeIndex[i+1] > endIndex: endIndex = self.nodeIndex[i+1]
                
                croppedX = self.x_res[self.nodeIndex[i]:self.nodeIndex[i+1]]
                croppedY = self.res[self.nodeIndex[i]:self.nodeIndex[i+1]]
                if self.graphics: plt.plot(croppedX, croppedY, '.-')
                if np.sum(croppedY) > 0:
                    sign = 1
                else:
                    sign = -1
                
                # Create guess, which should be the mid point of the zone
                guess = [1e25, (croppedX[0] + croppedX[-1])*0.5 ,0.001]
                fitG = bestGaus.select_Gaus_fit(croppedX,croppedY * sign, False, guess)
                fitResults = fitG.output()
                self.gauss_Fit_params.append(fitResults)
                xFit = np.arange(croppedX[0], croppedX[-1], 1e-5)
                yFit = self.gaus_selector(xFit, fitResults) * sign
#                self.newFit.extend( np.array(self.gaus_selector(croppedX, fitResults) * sign ))
                self.newX.extend(np.array(croppedX))
                if self.graphics: plt.plot(xFit,yFit)
                
                self.newFit += self.regionSelector_stepfuncs(self.nodesInRes[i], self.nodesInRes[i+1], self.regionCounter, fitResults, sign)
                
                
                #Need to create a function that creates cut offs for the curves
        FinalFormula = '  resFit = '
        for i in range(1,self.regionCounter+1):
            FinalFormula += 'finG{} + '.format(i)
        if self.displayToConsole: print FinalFormula[:-2]
        self.textOutput += FinalFormula[:-2]
        self.newFit = np.array(self.newFit)
#        print 'Fit params', self.gauss_Fit_params
#        print 'Start and end' , startIndex, endIndex
        if self.graphics:
            plt.xlim([self.x_res[startIndex], self.x_res[endIndex]])
            plt.show()        
    
            plt.plot(self.x_res, self.res, label = 'Input')
            plt.plot(self.x_res, self.newFit, lw = 3, label= 'new fit')
    #        plt.plot(self.newX, self.res[startIndex:endIndex] - self.newFit, label ='new residual')
            plt.legend()
            plt.show()
        return self.newFit, self.textOutput

    def printStepFnc(self, xc, k, up, number):
        if up: 
            if self.displayToConsole: print '  stpUp'+ str(number) +' = 0.5 + 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)
            self.textOutput += '  stpUp'+ str(number) +' = 0.5 + 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc) + '\n'
        else:
            if self.displayToConsole: print '  stpDwn'+ str(number) +' = 0.5 - 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)    
            self.textOutput += '  stpDwn'+ str(number) +' = 0.5 - 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc) + '\n'
        
    def regionSelector_stepfuncs(self, start, finish, number, gaussFunction, sign):
        if self.displayToConsole: print '#  Region Zone ', start, finish
        self.textOutput += '#  Region Zone ' + str(start) +' ' + str(finish) +'\n'
        
        #Print a line which is selecting each region
        yRegion = self.stepUp(self.x_res, start, self.steepness) * (1 - self.stepUp(self.x_res, finish, self.steepness)) * self.gaus_selector(self.x_res, gaussFunction) *sign
        if self.graphics: plt.plot(self.x_res, yRegion)
        
        self.printStepFnc(start, self.steepness, True, number)
        self.printStepFnc(finish, self.steepness, False,number)
        self.printGaussFormula_for_epoch(gaussFunction, number, sign)
        if self.displayToConsole: print '  finG{} = '.format(number) + 'stpUp{}'.format(number) + ' * stpDwn{}'.format(number) +' * rg{}'.format(number) 
        if self.displayToConsole: print 
        self.textOutput += '  finG{} = '.format(number) + 'stpUp{}'.format(number) + ' * stpDwn{}'.format(number) +' * rg{}'.format(number) +'\n\n'
        return yRegion
        
    def printGaussFormula_for_epoch(self, popt1, num, sign):
        outStr = ''
        outStr +='  rg{} = '.format(num) + '{0:.2e}'.format(popt1[1][0] * sign)
        outStr += ' * supergauss(x , {0:.4e}'.format (popt1[1][1]) 
        outStr += ', {0:.4e}'.format(popt1[1][2]**2) + ', '+ str(int(popt1[0])) +')'
        if self.displayToConsole: print outStr
        self.textOutput +=  outStr + '\n'
        
    def stepUp(self, x, xc = 0,  k = 1):
        return 0.5 + 0.5 * np.tanh(k * (x - xc))
        
    def gaus_selector(self, x,  inGausData):
        order = inGausData[0]
#        print 'Selecting Gaus of order: ' , order
        
        if order == 2:
#            print '2' , inGausData[1]
            return gaussian2(x, *inGausData[1])
        elif order == 4:
#            print '4', inGausData[1]
            return gaussian4(x, *inGausData[1])

        elif order == 6:
#            print '6'    , inGausData[1]
            return gaussian6(x, *inGausData[1])

        elif order == 8:
#            print '8'     , inGausData[1]  
            return gaussian8(x, *inGausData[1])

        

#    def createFittedFunction(self):
#        #Create the fitted function of gaussians and step functions
#        self.yoverall = self.gaus_selector(self.fitPopt[0]) 
##        self.yoverall += self.gaus_selector(self.fitPopt[1])

if __name__ == "__main__":
    resIn = np.loadtxt('testRes.txt')
    shock = 0.00015
    print len(resIn[:,0])
    
    resSolver = resFitter(resIn[:,0], resIn[:,1], shock, True, False)
    ResidualFit, Textoutput = resSolver.sections()
    print Textoutput

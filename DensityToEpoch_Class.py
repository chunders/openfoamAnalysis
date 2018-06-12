#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
       _ 
      /  |     | __  _ __  _
     /   |    /  |_||_|| ||
    /    |   /   |  |\ | ||_
   /____ |__/\ . |  | \|_|\_|
   __________________________ .
   
 ██████╗ ██████╗ ███████╗███╗   ██╗███████╗ ██████╗  █████╗ ███╗   ███╗
██╔═══██╗██╔══██╗██╔════╝████╗  ██║██╔════╝██╔═══██╗██╔══██╗████╗ ████║
██║   ██║██████╔╝█████╗  ██╔██╗ ██║█████╗  ██║   ██║███████║██╔████╔██║
██║   ██║██╔═══╝ ██╔══╝  ██║╚██╗██║██╔══╝  ██║   ██║██╔══██║██║╚██╔╝██║
╚██████╔╝██║     ███████╗██║ ╚████║██║     ╚██████╔╝██║  ██║██║ ╚═╝ ██║
 ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚═══╝╚═╝      ╚═════╝ ╚═╝  ╚═╝╚═╝     ╚═╝
                                                                       
████████╗ ██████╗                                                      
╚══██╔══╝██╔═══██╗                                                     
   ██║   ██║   ██║                                                     
   ██║   ██║   ██║                                                     
   ██║   ╚██████╔╝                                                     
   ╚═╝    ╚═════╝                                                      
                                                                       
███████╗██████╗  ██████╗  ██████╗██╗  ██╗                              
██╔════╝██╔══██╗██╔═══██╗██╔════╝██║  ██║                              
█████╗  ██████╔╝██║   ██║██║     ███████║                              
██╔══╝  ██╔═══╝ ██║   ██║██║     ██╔══██║                              
███████╗██║     ╚██████╔╝╚██████╗██║  ██║                              
╚══════╝╚═╝      ╚═════╝  ╚═════╝╚═╝  ╚═╝                              
                                                  
Created on Thu May 31 09:54:46 2018

@author: chrisunderwood

A Class for the fitting of the openfoam Density profile, so that the input file
can easily be altered.
"""
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('talk')    #paper, notebook, talk, poster
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec

import select_best_gaussFit_class as bestGaus


mm = True


def nearposn(array,value):
	#Find array position of value
    posn = (abs(array-value)).argmin()
    return posn

def gaussian(x, *params):
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

def stepUp(x, xc = 0,  k = 1):
    return 0.5 + 0.5 * np.tanh(k * (x - xc))

def printGaussFormula_for_epoch(popt1, num):
    print '  xg{} = '.format(num) + '{0:.2e}'.format(popt1[0]) + ' * supergauss(x , {0:.4e}'.format (popt1[1]) + ', {0:.4e}'.format(popt1[2] ** 2) + ', 6)'

def printStepFnc(xc, k):
    print 
    print '  stpUp = 0.5 + 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)
    print '  stpDwn = 0.5 - 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)


class Lineout_to_epoch():
    def __init__(self, DirectoryPath, fileName, invert, translation):
        self.shapeToFit = np.loadtxt(DirectoryPath + fileName)
        self.x = self.shapeToFit[:,0]*1e-3        #Convert to m, for input into epoch
        self.y = self.shapeToFit[:,1]
        if invert:
            self.y = self.y[::-1]
                
        self.ygrad = np.gradient(self.y)
        self.translate(translation)
        
        #Stepness of stepfunctions 
        self.stepness = 1e6   

        
    def translate(self, translation):
        #==============================================================================
        # Center on X = 0
        #==============================================================================
        self.shiftInX =  translation[0]
        self.x = self.x - self.shiftInX
        
        self.xScale = translation[1]
        self.x = self.x/self.xScale

    def plotInitProfile(self):
        plt.figure(figsize = (5,3.5))
        plt.plot(self.x,self.y)
        plt.vlines(0, 0, self.y.max(), 'r', linestyle = 'dashed')
        for i in self.x[self.indexs]:
            plt.vlines(i, 0, self.y.max(), 'g', linestyle = 'dashed')
        plt.title('Profile to plot')
        plt.show()
        
    def find_start_fin_shock(self, MoveShockToBottomMark):
        #This is either min or max depending on direction
        self.peakGrad = self.ygrad.min()
        
        self.shockFront =  self.x[nearposn(self.ygrad, self.peakGrad)]
        
#        print 
#        print 'Peak gradient: ', self.ygrad.max(), nearposn(self.ygrad, self.ygrad.max())
#        print 'Shock pos: ', self.shockFront
        
        yabs = abs(self.ygrad)
        start = np.argmax(yabs >  abs(self.peakGrad*0.05))
        finish = np.argmax(yabs[::-1] >  abs(self.peakGrad*0.1))
    #    print x[start], x[::-1][finish]
        self.indexs = [start, nearposn(self.ygrad, self.peakGrad), len(self.x) - finish]
        self.indexs[1] = self.indexs[1] + MoveShockToBottomMark
#        self.plotInitProfile()
        self.stepPos = self.x[self.indexs[1]]

    
    def printIndexes(self):
        print 'Index of start, shock and end: ', self.indexs

    def fit_main_superGauss(self, guess):
        self.fitPopt = []
        for i in range(len(self.indexs) - 1):
            #==========================================================================
            #     Fitting the two main gaussians
            #==========================================================================
            xcrop = self.x[self.indexs[i]: self.indexs[i+1]]
            ycrop = self.y[self.indexs[i]: self.indexs[i+1]]
            
            guess[i][1] = (xcrop[0] + xcrop[-1]) * 0.5
            
            fitG = bestGaus.select_Gaus_fit(xcrop, ycrop, False, guess[i])
            print 'Output of new fitting function'
            a, b =  fitG.output()
            print a, b
            popt, pcov = curve_fit(gaussian6, xcrop,ycrop , p0=guess[i]) #, bounds=(0, [1e27, 20, 20]))
            if popt[2] < 0:
                print 'Width is negative, taking abs value'
                popt[2] = abs(popt[2])
                poptInit = popt
                popt, pcov = curve_fit(gaussian6, xcrop,ycrop , p0=popt)
                print 'ratio of fits', poptInit / popt
                
        #    plt.plot(xcrop,ycrop)   
        #    plt.plot(xcrop, gaussian6(xcrop, *popt), '--r', lw = 3)
        #    plt.plot(x, gauss6(x,popt[0], popt[1], popt[2]), '-.b')    
        #    plt.vlines(popt[1], 0, 1e26)
        #    fact = 2
        #    plt.vlines(popt[1] + popt[2]**(fact), 0, 1e26)
        #    plt.vlines(popt[1] - popt[2]**(fact), 0, 1e26)
        #    plt.show()
                
#            print 'Fit:'
#            for p, e in zip( popt, np.sqrt(np.diag(pcov))/popt):
#                print '\t' + str(p) + ' +/- ' + str(e)
#            
        
            self.fitPopt.append(popt)
    
        self.fitPopt = np.array(self.fitPopt)     
        
    def createFittedFunction(self):
        #Create the fitted function of gaussians and step functions
        self.yoverall = gaussian6(self.x, *self.fitPopt[0])  * (1-stepUp(self.x, xc = self.stepPos, k = self.stepness)) 
        self.yoverall += gaussian6(self.x, *self.fitPopt[1]) * (stepUp(self.x, self.stepPos, k = self.stepness))
        
        self.res = (self.y - self.yoverall) #/self.yoverall
        
        self.res = self.res[self.indexs[0]:self.indexs[2]]
        self.x_res = self.x [self.indexs[0]:self.indexs[2]]
        
    def plotFit_vs_input(self):
        #==============================================================================
        # plotting the resulting function 
        #==============================================================================
        plt.figure(figsize = (8,8))
    
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        
        if mm:
            self.x = self.x *1e3
            self.x_res = self.x_res *1e3
        
        ax1 = plt.subplot(gs[0])
    #    ax1.plot(xFit, yfit2, label = 'Attempted model')
        ax1.plot(self.x, self.yoverall, 'b', lw = 3, label = 'Fit')
        ax1.plot(self.x,self.y, 'tab:orange', label = 'data to fit: openfoam')
        ax1.set_xlim([self.x[0], self.x[-1]])
        ax1.legend()
        ax1.vlines(0, 0, self.y.max(), colors='r', linestyle='dashed')
        
        ax1.yaxis.grid(color='gray', linestyle='dashed')
        ax1.xaxis.grid(color='gray', linestyle='dashed')
        
        for i in self.x[self.indexs]:
            ax1.vlines(i, 0, self.y.max(), 'g', linestyle = 'dashed')
        
        ax2 = plt.subplot(gs[1], sharex = ax1)
        
            
        ax2.plot(self.x, self.yoverall, 'b',  lw = 2)
        ax2.set_xlim([self.x[0], self.x[-1]])
        
        ax2.set_ylabel('Density', color='b')
        ax2.tick_params('y', colors='b')
        ax3 = ax2.twinx()
        
        
        ax3.set_axisbelow(True)
        ax3.yaxis.grid(color='red', linestyle='dashed')
        ax3.xaxis.grid(color='red', linestyle='dashed')
        ax3.set_xlim([self.x[0], self.x[-1]])
        ax3.plot(self.x_res, self.res,  'r') 
        ax3.set_ylabel('Residual', color='r')
        ax3.tick_params('y', colors='r')
        if mm:
            ax2.set_xlabel('Distance (mm)')
            self.x = self.x *1e-3
            self.x_res = self.x_res *1e-3
        else:
            ax2.set_xlabel('Distance (m)')
        
    
    #    plt.savefig('FittingOpenFoamToEPOCH.png', dpi=300)
#        plt.show()
        
    

        
    def displayTextForEpoch(self, name):
        yEnd = len(self.yoverall) - next((i for i, val in enumerate((self.yoverall[::-1] - self.yoverall[::-1].min() / self.yoverall[::-1].max())) if val > (self.yoverall.max() - self.yoverall.min()) / 40), None)
        print; print 'Distance at end -> ', yEnd, 'is : ', self.x[yEnd], 'm or ' , self.x[yEnd] * 1e6, 'um'
        print 'Shock index     -> ', self.indexs[1], 'is : ', self.x[self.indexs[1]], 'm or ', self.x[self.indexs[1]] * 1e6, 'um'
        
        print; print 
        print 'In constant block'; print
        print '  # ' + name
        for popt, i in zip(self.fitPopt, range(1, len(self.fitPopt)+1)):
            printGaussFormula_for_epoch(popt, i)
        printStepFnc(self.stepPos, self.stepness)
        print '  xDens = xg1 * stpDwn + xg2 * stpUp'
        print; print 'in control block'
        time = (self.x[yEnd]) / 3e8
        print 
        print '  t_end = {:.3e}'.format(time*1.01)
                
def run(folderPath, fileName, invert, MoveShockToBottomMark,
        shiftInX, xScale, guess, name):
    translate = [0.0, xScale]
    
    #Locate the shock front and them shift to shiftInX infront of it
    dp = Lineout_to_epoch(folderPath, fileName, invert, translate)
    dp.find_start_fin_shock(MoveShockToBottomMark)
    
    #Moving the whole data to be just in front of shock at zero
    translate = [dp.shockFront-shiftInX / xScale, 1]
    dp.translate(translate)
    dp.find_start_fin_shock(MoveShockToBottomMark)

#    dp.printIndexes()
    dp.fit_main_superGauss(guess)
    dp.createFittedFunction()
    dp.plotFit_vs_input()
    dp.displayTextForEpoch(name)
    plt.savefig('densityMapping__' + folderPath.split('/')[-2] + '__' + fileName.split('.')[0] + '.png')
    plt.show()


if __name__ == "__main__": 
    folderPath = '/Volumes/CIDU_passport/openFOAM/Line Out Data/Blade_lineout/h80_0mm/'
    fileName = '45lineout.txt'
    invert = True
    MoveShockToBottomMark = 15
    
    
#    folderPath = '/Volumes/CIDU_passport/openFOAM/Line Out Data/1mm_Above_Blade_parallel/'
#    folderPath += 'h50_t-50/'
#    fileName = '95lineout.txt'
#    invert = False
#    MoveShockToBottomMark = 15

    shiftInX =  1.5e-3
    xScale = 10
    guess = [[2e26, 1e-3  /xScale, 10e-3],
             [1e26, 7.5e-3/xScale, 10e-3]]

    run(folderPath, fileName, invert, MoveShockToBottomMark,
        shiftInX, xScale, guess, '')
    


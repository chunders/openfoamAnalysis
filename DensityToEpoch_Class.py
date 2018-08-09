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
import matplotlib.gridspec as gridspec

import select_best_gaussFit_class as bestGaus
if False:
    import foamToEpoch_ResidualFit_Class_working as resFit
else:
    import foamToEpoch_ResidualFit_Class as resFit



SNS = False
if SNS:
    import seaborn as sns
    sns.set_context('poster')    #paper, notebook, talk, poster

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

def gaussian8(x, *params):
    #Gaussian function
    A = params[0]
    x0 = params[1]
    c = params[2]
    return A*np.exp(-((x-x0)/(c**2))**8)

def supergauss(x, x0, c, order):
    return np.exp(-((x-x0)/(c**2))**order)

def stepUp(x, xc = 0,  k = 1):
    return 0.5 + 0.5 * np.tanh(k * (x - xc))

def printGaussFormula_for_epoch(popt1, num):
    outStr = ''
    outStr +='  xg{} = '.format(num) + '{0:.2e}'.format(popt1[1][0])
    outStr += ' * supergauss(x , {0:.4e}'.format (popt1[1][1]) 
    outStr += ', {0:.4e}'.format(popt1[1][2]**2) + ', '+ str(int(popt1[0])) +')'
    print outStr

def printStepFnc(xc, k):
    print 
    print '  stpUp = 0.5 + 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)
    print '  stpDwn = 0.5 - 0.5*tanh({0:.1e}'.format(k) + ' * (x - {0:.4e}))'.format(xc)


class Lineout_to_epoch():
    def __init__(self, DirectoryPath, fileName, invert, translation):
        self.shapeToFit = np.loadtxt(DirectoryPath + fileName)
        self.fileName = fileName
        self.folderPath = DirectoryPath
        self.x = self.shapeToFit[:,0]*1e-3        #Convert to m, for input into epoch
        self.y = self.shapeToFit[:,1]
        if invert:
            self.y = self.y[::-1]
                
        self.ygrad = np.gradient(self.y)
        self.translate(translation)
        self.xScale = translation[1]
        #steepness of stepfunctions 
        self.steepness = 1e6   

        
    def translate(self, translation):
        #==============================================================================
        # Center on X = 0
        #==============================================================================
        self.shiftInX =  translation[0]
        self.x = self.x - self.shiftInX
        
        self.xScale2 = translation[1]
        self.x = self.x/self.xScale2

    def plotInitProfile(self):
        plt.figure(figsize = (5,3.5))
        plt.plot(self.x,self.y)
        plt.vlines(0, 0, self.y.max(), 'r', linestyle = 'dashed')
        for i, col in zip(self.x[self.indexs], ['g', 'k', 'b']):
            plt.vlines(i, 0, self.y.max(), col, linestyle = 'dashed')
        plt.title('Profile to plot')
        plt.show()
        
    def find_start_fin_shock(self, MoveShockToBottomMark):
        #This is either min or max depending on direction
        print 'finding indexes'
        self.peakGrad = self.ygrad.min()
        
        self.shockFront =  self.x[nearposn(self.ygrad, self.peakGrad)]
        shockIndex = nearposn(self.ygrad, self.peakGrad)
#        print 
#        print 'Peak gradient: ', self.ygrad.max(), nearposn(self.ygrad, self.ygrad.max())
#        print 'Shock pos: ', self.shockFront
        
        yabs = abs(self.ygrad)
        start = np.argmax(yabs >  abs(self.peakGrad*0.05))
        finish = np.argmax(yabs[::-1] >  abs(self.peakGrad*0.1))
#        print finish, len(yabs)
        self.indexs = [start, shockIndex, len(self.x) - finish]
        if finish == 0:
#            print 'Last index Error'
            finish = shockIndex + nearposn(self.y[shockIndex:], self.y[shockIndex:].min())
            print nearposn(self.y[shockIndex:], self.y[shockIndex:].min())
#            plt.plot(self.x[shockIndex:], self.y[shockIndex:])
#            plt.vlines(self.x[finish], 0, self.y.max())
            self.indexs = [start, shockIndex, finish]

#        print finish, len(yabs)

    #    print x[start], x[::-1][finish]
        self.indexs[1] = self.indexs[1] + MoveShockToBottomMark
        self.plotInitProfile()
        self.stepPos = self.x[self.indexs[1]]
        print self.indexs


    
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
            sgOrder, poptBestGaus =  fitG.output()
            print sgOrder, poptBestGaus
            
#            popt, pcov = curve_fit(gaussian6, xcrop,ycrop , p0=guess[i]) #, bounds=(0, [1e27, 20, 20]))
#            if popt[2] < 0:
#                print 'Width is negative, taking abs value'
#                popt[2] = abs(popt[2])
#                poptInit = popt
#                popt, pcov = curve_fit(gaussian6, xcrop,ycrop , p0=popt)
#                print 'ratio of fits', poptInit / popt
                
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
        
            self.fitPopt.append([sgOrder, poptBestGaus])
    
        print 'Output of fit'
        print self.fitPopt
#        self.fitPopt = np.array(self.fitPopt)     
        
    def gaus_selector(self, inGausData):
        order = inGausData[0]
        print 'Selecting Gaus of order: ' , order
        
        if order == 2:
            print '2' , inGausData[1]
            return gaussian(self.x, *inGausData[1])
        elif order == 4:
            print '4', inGausData[1]
            return gaussian4(self.x, *inGausData[1])

        elif order == 6:
            print '6'    , inGausData[1]
            return gaussian6(self.x, *inGausData[1])

        elif order == 8:
            print '8'     , inGausData[1]  
            return gaussian8(self.x, *inGausData[1])

        

    def createFittedFunction(self):
        #Create the fitted function of gaussians and step functions
        self.yoverall = self.gaus_selector(self.fitPopt[0])  * (1-stepUp(self.x, xc = self.stepPos, k = self.steepness)) 
        self.yoverall += self.gaus_selector(self.fitPopt[1])  * (stepUp(self.x, self.stepPos, k = self.steepness))
        
        self.res = (self.y - self.yoverall) #/self.yoverall
        
        self.res = self.res[self.indexs[0]:self.indexs[2]]
        self.x_res = self.x [self.indexs[0]:self.indexs[2]]
        
    def plotFit_vs_input(self):
        #==============================================================================
        # plotting the resulting function 
        #==============================================================================
        plt.figure(figsize = (8,8))
        if SNS: sns.set_style('white')
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
        
        ax1.yaxis.grid(color='gray', linestyle='dashed' , lw = 0.5)
        ax1.xaxis.grid(color='gray', linestyle='dashed', lw = 0.5)
        
        #   Add lines to show where the main points are
        for i in self.x[self.indexs]:
            ax1.vlines(i, 0, self.y.max(), 'g', linestyle = 'dashed')
        
        #   Create axis with the scale
        ax2 = plt.subplot(gs[1], sharex = ax1)
        #   Plot the function
        ax2.plot(self.x, self.yoverall, 'b',  lw = 2)
        ax2.set_xlim([self.x[0], self.x[-1]])
        
        ax2.set_ylabel('Density', color='b')
        ax2.tick_params('y', colors='b')
        
        #   Twin axis to put residual on same plot
        ax3 = ax2.twinx()      
        
        ax3.set_axisbelow(True)
        ax3.yaxis.grid(color='red', linestyle='dashed')
        ax3.xaxis.grid(color='red', linestyle='dashed')
        ax3.set_xlim([self.x[0], self.x[-1]])
        ax3.plot(self.x_res, self.res,  'r', lw = 1.5) 
        ax3.plot(self.x[self.indexs[0]:self.indexs[2]], 
                 self.y[self.indexs[0]:self.indexs[2]] - self.yoverall[self.indexs[0]:self.indexs[2]],
                 'g', lw = 0.9) 
        ax3.plot(self.nodes * 1000, np.zeros(len(self.nodes)), '.')

        ax3.set_ylabel('Residual', color='r')
        ax3.tick_params('y', colors='r')
        if mm:
            ax2.set_xlabel('Distance (mm)')
            self.x = self.x *1e-3
            self.x_res = self.x_res *1e-3
        else:
            ax2.set_xlabel('Distance (m)')
            
        #   Save the residual to work with later
        np.savetxt(self.folderPath + 'residualToFit_' +self.fileName + '.txt', 
                   np.c_[self.x_res, self.res])
        plt.suptitle('OpenFoam lineout to Epoch input deck function')
#        plt.tight_layout()
    #    plt.savefig('FittingOpenFoamToEPOCH.png', dpi=300)
#        plt.show()
        
    

        
    def displayTextForEpoch(self, name, resFitText):
        yEnd = len(self.yoverall) - next((i for i, val in enumerate((self.yoverall[::-1] - self.yoverall[::-1].min() / self.yoverall[::-1].max())) if val > (self.yoverall.max() - self.yoverall.min()) / 40), None)
        print; print 'Distance at end -> ', yEnd, 'is: ', self.x[yEnd], 'm or ' , self.x[yEnd] * 1e6, 'um'
        print 'Shock index     -> ', self.indexs[1], 'is: ', self.x[self.indexs[1]], 'm or ', self.x[self.indexs[1]] * 1e6, 'um'
        
        print; print 'Fit params to input'
        for i in self.fitPopt:
            print i
        
        print; print 
        print 'In constant block'; print
        print '  # ' + name + 'Scaling: ' + str(self.xScale)
        for popt, i in zip(self.fitPopt, range(1, len(self.fitPopt)+1)):
            printGaussFormula_for_epoch(popt, i)
            
        printStepFnc(self.stepPos, self.steepness)
        print '  xDens = xg1 * stpDwn + xg2 * stpUp'
        print 
        print resFitText
        print 
        print '  xDensProfile = xDens + resFit'
        print 
        
        print; print 'in control block'
        time = (self.x[yEnd]) / 3e8
        print 
        print '  t_end = {:.3e}'.format(time*1.01)
        
#==============================================================================
#     This function should become redundant with the new file being created
#     to look into this
#==============================================================================
    def inspect_residual(self, regExtender, FitBeyond):
        
        extraFit = resFit.resFitter(self.x_res, self.res, self.shockFront, False, False, regExtender, FitBeyond)
        newFit, resFitText, self.nodes = extraFit.sections()

        np.savetxt('testRes.txt', np.c_[self.x_res, self.res])

        self.yoverall[self.indexs[0]:self.indexs[-1]] += newFit
        return resFitText
        
#        plt.plot(self.x_res, self.res)
#        crossings = []
#        for i in range(len(self.res)-1):
#            if self.res[i] * self.res[i+1] < 0:
#                crossings.append(i)
#        
#        self.nodesInRes = self.x_res[crossings] 
#        plt.plot(self.nodesInRes, np.zeros(len(crossings)), 's')
#        
#        plt.show()
#        print self.nodesInRes
#        print self.shockFront
        
                
def run(folderPath, fileName, invert, MoveShockToBottomMark,
        shiftInX, xScale, guess, name, regExtender = 0, FitBeyond = 5e-05):
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
    resFitText = dp.inspect_residual(regExtender, FitBeyond)

    dp.plotFit_vs_input()
    plt.savefig(folderPath + 'densityMapping__' + folderPath.split('/')[-2] + '__' + fileName.split('.')[0] + '.png',
                dpi = 300)
    plt.show()
    
    dp.displayTextForEpoch(name, resFitText)



if __name__ == "__main__": 
    folderPath = '/Volumes/CIDU_passport/openFOAM/Line Out Data/Blade_lineout/h80_0mm/'
    folderPath = '/Volumes/GoogleDrive/My Drive/HydroSimulations_openFOAM/Line Out Data/Blade_lineout/h80_0mm/'
    fileName = '45lineout.txt'
#    fileName = '65lineout.txt'
    invert = True
    MoveShockToBottomMark = 15
    
    
#    folderPath = '/Volumes/CIDU_passport/openFOAM/Line Out Data/1mm_Above_Blade_parallel/'
#    folderPath = '/Users/chrisunderwood/Downloads/1mm_Above_Blade_parallel/'

#    folderPath += 'h50_t-65/'
#    fileName = '95lineout.txt'
#    fileName = '96lineout.txt'
#    invert = False
#    MoveShockToBottomMark = 15

    shiftInX =  1.5e-3
    xScale = 10
    guess = [[2e26, 1e-3  /xScale, 10e-3],
             [1e26, 7.5e-3/xScale, 10e-3]]

    run(folderPath, fileName, invert, MoveShockToBottomMark,
        shiftInX, xScale, guess, folderPath.split('/')[-2])
    
    #   if needing to check things
    x = np.arange(-0.001, 0.002, 0.000001)
    


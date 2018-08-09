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
                                                                       
Created on Fri Jun  1 10:27:28 2018

@author: chrisunderwood

Average Lineout from Shock Analysis.

Create plots looking at the different position of the blade 
and how they effect the results

!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!

This file needs the following file to have been ran FIRST

/Users/chrisunderwood/Documents/Gas Jet Simulation - OpenFOAM/
Lineout_Density_profile_analysis_Class.py


"""
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import gridspec
import math
import seaborn as sns

sns.set_context('poster')
figSize = (10,10)

#==============================================================================
# Read in the list of subfolders that contain all the data that is to be 
# analysed.
# The file read in may have to be altered to look at different file outputs
#==============================================================================
class AveLineout():
    def __init__(self, Subfolders):
        self.subFolders = Subfolders
        
        #Find the name of the subfolder from the list
        self.idFolder = []
        for sub in self.subFolders:
            arr = sub.split('/')
            self.idFolder.append(arr[-2])          
#        print 'Folder id: ', self.idFolder
    
        # Find the positon of the blade from the folder name
        self.bladePos = []
        for name in self.idFolder:
            namesplit = name.split('_')
            pos = np.zeros(2)
            for n in namesplit:
                if n.startswith('t'):
                    pos[0] = int(n[1:])
                if n.startswith('h'):
                    pos[1] = int(n[1:])
            self.bladePos.append(pos)
        self.bladePos = np.array(self.bladePos)
                    
#        print 'Blade Positions'
#        print self.bladePos
        
        # Read in the data file for each subfolder
        self.d = []
        self.regionDims = []        
        for sub , idF in zip(self.subFolders, self.idFolder):
#            print sub
            fileName = sub + 'Average_lineout_of_shock_' + idF + '.txt'
#            print fileName
            self.d.append( np.loadtxt(fileName))
            self.regionDims.append([0,0,0,0])
        self.regionDims = np.array(self.regionDims, dtype = 'float')
        #Scan through and work out different heights
        self.height = self.sortPosition(1)
        #Scan through and work out different heights
        self.trans = self.sortPosition(0)        
            
#==============================================================================
# A function that plots all the read in data on to the same graph
#==============================================================================
    def plotAllLines(self, show = True):
        c=iter(plt.cm.nipy_spectral(np.linspace(0.2,0.9,len(self.d))))
        for data, lab in zip(self.d, self.idFolder):
            col = next(c)
            plt.plot(data[:,0], data[:,1], label = lab, color = col)
        plt.xlabel('Distance (mm)')
        plt.ylabel(r'Number Density $(m^{-3})$')
        if show:
            plt.legend()
            plt.show()
            
#==============================================================================
# From the file name sorts them into there identifying number and which files
# are in that. The h_or_t represents which number to go from
#                 h_or_t = 1 --> Height
#                 h_or_t = 0 --> Translation
#==============================================================================
    def sortPosition(self, h_or_t):
        scanPos = []
        for i in range(len(self.bladePos)):
            knownscanPos = 0
            if not len(scanPos) == 0:
                for h in scanPos:
                    if self.bladePos[i][h_or_t] == h[0]:
                        knownscanPos = 1
                    
            if knownscanPos == 0:
                print 'Adding to list: ' , self.bladePos[i][h_or_t]
                # Append the scanPos and an array to fill with the indexes
                # of the data for that scanPos
                scanPos.append([self.bladePos[i][h_or_t],[]])
                
            for h in scanPos:
                if self.bladePos[i][h_or_t] == h[0]:
                    h[1].append(i)
        print 'The number of different options: ', len(scanPos)
        self.nos_of_params = len(scanPos)
        return scanPos
        
#==============================================================================
# Plot the data looking at the different translations but constant height
#==============================================================================
    def plot_constH(self):
        plt.figure(figsize = figSize) 
        gs = gridspec.GridSpec(nrows = int(math.ceil(self.nos_of_params / 2)), ncols = 2) 

        for plotHeight, gsPos in zip(self.height, gs):
            ax = plt.subplot(gsPos)
            for index in plotHeight[1]:
                ax.plot(self.d[index][:,0], self.d[index][:,1], 
                         label = self.bladePos[index][0]*0.1)
                
            ax.set_title('Height of blade: ' + str(plotHeight[0]/ 10 )+ ' mm')
            ax.legend(title = 'Translation (mm)')
        plt.tight_layout()
        plt.savefig(mainDir + 'Const_height_lineout_comparison.png', dpi = 250)
        plt.show()
          
#==============================================================================
# Plot the data looking at the different height but constant translation
#==============================================================================
    def plot_constT(self):       
        plt.figure(figsize = figSize) 
        gs = gridspec.GridSpec(nrows = int(math.ceil(self.nos_of_params / 2)), ncols = 2) 

        for plotTrans , gsPos in zip(self.trans, gs):
            ax = plt.subplot(gsPos)

            for index in plotTrans[1]:
                ax.plot(self.d[index][:,0], self.d[index][:,1], 
                         label = self.bladePos[index][1]*0.1)
                
            ax.set_title('Tranlation of blade: ' + str(plotTrans[0]/ 10 )+ ' mm')
            ax.legend(title = 'Height (mm)')
        
        plt.tight_layout()
        plt.savefig(mainDir + 'Const_trans_lineout_comparison.png', dpi = 250)
        plt.show()
        
    def print_index(self):
        for i in range(len(self.idFolder)):
            print i, self.idFolder[i]
    
    def get_data(self, index):
        print self.idFolder[index]
        return self.d[index]
        
    def get_trans(self):
        print 'Sorted by translation'
        for i, op in zip(range(len(self.trans)), self.trans):
            print i, op
        return self.trans
        
    def get_height(self):
        print 'Sorted by height'
        for i, op in zip(range(len(self.height)), self.height):
            print i, op
        return self.height

    def plotOneOption(self, option, trans0height1):
        print 'Name ', option[0]
        
        c=iter(plt.cm.tab10(np.linspace(0.0,0.9,len(option[1]))))
        
        for index in option[1]:
            col = next(c)
            plt.plot(self.d[index][:,0], self.d[index][:,1], 
                         label = self.bladePos[index][trans0height1]*0.1,
                         color = col)
        plt.xlabel('Distance (mm)')
        plt.ylabel(r'Number Density $(m^{-3})$')
        if trans0height1 == 0: 
            title = 'Height ' 
        else:
            title = 'Trans '
        plt.title(title + str(option[0]*0.1) + 'mm')
        if trans0height1 == 0: 
            title = 'Trans '
        else:
            title = 'Height '
        plt.legend(title = title + '(mm)')
        
    def findRegions(self, index):
        x = self.d[index][:,0]        
        line = self.d[index][:,1]
        grad = np.gradient(line)
#        N = 5; 
#        filt = np.convolve(grad , np.ones((N,))/N, mode='valid')
        fig, ax = plt.subplots(1,1,figsize=(5,5))

        fitProf = Lineout_Pos_Locator(x, line)
        pos = fitProf.find_start_fin_shock(25)
        print 'Indexes', pos
        print 'Xpos: ', x[pos]
        shapeMarkers = x[pos]
        
        highRegion = shapeMarkers[1] - shapeMarkers[0]
        lowRegion = shapeMarkers[3]-shapeMarkers[1]
        
        
        print 'L1:', highRegion
        print 'L2: ', lowRegion
        errLow = (shapeMarkers[4]-shapeMarkers[1]) - lowRegion
        errHigh = (shapeMarkers[1] - shapeMarkers[5]) - highRegion
        print 'L1err: ',  errHigh
        print 'L2err: ', errLow
        
        
        
        self.regionDims[index] = np.array([shapeMarkers[1] - shapeMarkers[0],shapeMarkers[3]-shapeMarkers[1] , errLow, errHigh])
        plt.plot(x, line, 'b', lw = 3)
        plt.ylabel('Density')
        for index in pos[:2], pos[3] :
            plt.vlines(x[index], 0, 1.5e26, colors='green')
            
        ax2 = ax.twinx()
        ax2.plot(x, grad, lw = 0.75)
        ax2.plot(x, abs(grad), lw = 0.75)
        '''
            It may be possible to do some kind of turning point
            fit to work out the end of the flat region
        '''         
        plt.ylabel('Density Grad')
        print 'Y vals at boundaries: ', line[pos]
#        ax2.plot(x[:len(filt)],filt)
        plt.show()
        
        
        
        
        
def nearposn(array,value):
	#Find array position of value
    posn = (abs(array-value)).argmin()
    return posn

class Lineout_Pos_Locator():
    def __init__(self, X, Y):
        self.x = X        #Convert to m, for input into epoch
        self.y = Y              
        self.ygrad = np.gradient(self.y)
        
    def find_start_fin_shock(self, MoveShockToBottomMark):
        #This is either min or max depending on direction
        self.peakGrad = self.ygrad.min()
        self.shockIndex =  nearposn(self.ygrad, self.peakGrad)
        self.shockFront =  self.x[self.shockIndex]
        
        print 
        print 'Peak gradient: ', self.ygrad.max(), nearposn(self.ygrad, self.ygrad.max())
        print 'Shock pos: ', self.shockFront
        
        yabs = abs(self.ygrad)
        start = np.argmax(yabs >  abs(self.peakGrad*0.2))
        finish = np.argmax(yabs[::-1] >  abs(self.peakGrad*0.1))
        
        startErr = np.argmax(yabs >  abs(self.peakGrad*0.25))

        
    #    print x[start], x[::-1][finish]
        self.indexs = [start, nearposn(self.ygrad, self.peakGrad), len(self.x) - finish]
        self.indexs[1] = self.indexs[1] + MoveShockToBottomMark
        
        self.ycrop = self.y[self.indexs[1]:]
        tenPercentDrop = np.argmax(self.ycrop[0]*0.9 >  abs(self.ycrop))
        twentyPercentDrop = np.argmax(self.ycrop[0]*0.85 >  abs(self.ycrop))

        print 'NewEnd pos: ', tenPercentDrop
        print 'or ', self.x[self.indexs[1]+tenPercentDrop]
        
        self.stepPos = self.x[self.indexs[1]]
        self.indexs.append(self.indexs[1]+tenPercentDrop)
        self.indexs.append(self.indexs[1]+twentyPercentDrop)
        self.indexs.append(startErr)
        return self.indexs


if __name__ == "__main__":
    mainDir = '/Volumes/CIDU_passport/openFOAM/Line Out Data/1mm_Above_Blade_parallel/'
#    mainDir = '/Users/chrisunderwood/Downloads/1mm_Above_Blade_parallel/'

    listSubFolders =  [x[0]+ '/' for x in os.walk(mainDir)][1:]
    print listSubFolders
    
    #lines = AveLineout(listSubFolders[:8])
    lines = AveLineout(listSubFolders)
    
    #lines.plotAllLines()
    lines.plot_constH()
    #print ; print
    #lines.plot_constT()
    
    print;print 
    trans = lines.get_trans()
    height = lines.get_height()
    
    lines.plotOneOption([50.0, [0, 1,2, 3,  5,  7]], 0)
#    lines.plotOneOption(height[0], 0)

    plt.tight_layout()
    plt.savefig(mainDir + 'Altering_Trans.png')
    plt.show()
    
    lines.plotOneOption(trans[7], 1)
    plt.tight_layout()
    plt.savefig(mainDir + 'Altering_Height.png')
    plt.show()
    
    for i in range(8):
        lines.findRegions(i)
        
    plotFit = lines.regionDims[height[0][1]]
    def removeZeros2d(arr):
        out = []
        for a in arr:
            print a[0]
            if a[0] > 0.0:
                out.append(a)
        return np.array(out)
    #plotFit = removeZeros2d(plotFit)
    
    #plt.plot(plotFit[:,0], plotFit[:,1], 'o-')
    #plt.xlabel('Initial Region Length')
    #plt.ylabel('Accelerating Region Length')
    #
    #for label, x, y in zip(lines.bladePos[height[0][1]] , plotFit[:, 0], plotFit[:, 1]):
    #    plt.annotate(
    #        str(label[0]) + ',' + str(label[1]),
    #        xy=(x, y), xytext=(30, 20),
    #        textcoords='offset points', ha='right', va='bottom',
    #        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
    #        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))
    #plt.show()
    
        
    capSize = 5
    capWidth = 2
    plt.figure(figsize=(8,8))
    ax = plt.subplot()
    sns.set_style('whitegrid')
    (_, caps, _) = ax.errorbar(lines.bladePos[height[0][1]][:,0] * 0.1,
                                               plotFit[:, 1] / plotFit[:, 0],
                  yerr = (plotFit[:, 1] + plotFit[:, 2])*0.5 / plotFit[:,0],
                  capsize=capSize)
    for cap in caps:
        cap.set_markeredgewidth(capWidth)
          
    
    plt.xlabel('Blade horizontal pos from centre of nozzle (mm)')
    plt.ylabel('Ratio of accelerating region to initial region')
    
    def lin(x, *params):
        return x* params[0] + params[1]
    
    from scipy.optimize import curve_fit
    popt, pcov = curve_fit(lin, lines.bladePos[height[0][1]][:,0] * 0.1,
                           plotFit[:, 1] / plotFit[:, 0] , p0=[-0.5, 0]) #, bounds=(0, [1e27, 20, 20]))
    perr = np.sqrt(np.diag(pcov))
    xfit = np.arange(-5.5, -0.9, 0.1)
    yfit = lin(xfit, *popt)
    ax.plot(xfit, yfit)
    ax.text(-4.5, 2.1, 
             'y = ({0:.2f}'.format(popt[0])+ '+/-{0:.2f}'.format(perr[0])+ ')*x + ({0:.2f}'.format(popt[1])+ '+/-({0:.2f}'.format(perr[1])+ ')',
             bbox=dict(facecolor='white', alpha=0.5))
    plt.show()
    
    
    
    #==============================================================================
    # A plot that shows the trends of both regions, high and accelerating
    #==============================================================================
    sns.set_style('white')
    plt.figure(figsize=(8,8))
    xfit = np.arange(-7, -0.9, 0.1)
    
    #High Density Region
    x = lines.bladePos[height[0][1]][:,0] * 0.1
    y = plotFit[:, 1]
    (_, caps, _) = plt.errorbar(x,y,
                  yerr = (plotFit[:, 1] + plotFit[:, 2])*0.5 ,
                  color= 'b',
                  lw = 0.75, 
                  fmt = 'o',
                  capsize=capSize,
                  label = 'Low Density Region')
    for cap in caps:
        cap.set_markeredgewidth(capWidth)
        
    popt, pcov = curve_fit(lin, x,y, p0=[0.5, 0])
    yfit = lin(xfit, *popt)
    plt.plot(xfit, yfit, 'b')
    plt.text(-6.8, 15, 
             'y = ({0:.2f}'.format(popt[0])+ '+/-{0:.2f}'.format(perr[0])+ ')*x + ({0:.2f}'.format(popt[1])+ '+/-({0:.2f}'.format(perr[1])+ ')',
             bbox=dict(facecolor='blue', alpha=0.5))
    
    #Low Density Region
    x = lines.bladePos[height[0][1]][:,0] * 0.1
    y = plotFit[:, 0]
    (_, caps, _) = plt.errorbar(x,y,
                  yerr = (plotFit[:, 0] + plotFit[:, 3])*0.5,
                  color= 'r',
                  capsize=capSize,
                  lw = 0.75, 
                  fmt = 'o',
                  label = 'High Density Region')
    for cap in caps:
        cap.set_markeredgewidth(capWidth)
    
    popt, pcov = curve_fit(lin, x,y, p0=[-0.5, 0])
    yfit = lin(xfit, *popt)
    plt.plot(xfit, yfit, 'r')
    plt.text(-6.8, 13.7, 
             'y = ({0:.2f}'.format(popt[0])+ '+/-{0:.2f}'.format(perr[0])+ ')*x + ({0:.2f}'.format(popt[1])+ '+/-({0:.2f}'.format(perr[1])+ ')',
             bbox=dict(facecolor='red', alpha=0.5))
             
             
    plt.legend()
    plt.xlabel('Blade Translation from Centre (mm)')
    plt.ylabel('Size of Region (mm)')


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
  
Created on Thu May 24 11:21:14 2018

@author: chrisunderwood

Create a Class python script for looking at the different lineouts with the 
blades

Plots the 2D evolution
Saves the average profile during this time from the times step start point

Has the option of producing line outs at individual time steps too.

This file has to be run before
/Users/chrisunderwood/Documents/Gas Jet Simulation - OpenFOAM/
Averaged_Lineout_Data_Analysis.py


"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib import gridspec

class LineoutAnalysis():
    def __init__(self, DirectoryPath, ma = False):
        self.folderPath = DirectoryPath

        self.csv_list, self.timeSteps = self.FilesInFolder()
        self.mach = ma
        self.assignSpreadsheetHeader()
        self.identifingFolder =  self.folderPath.split('/')[-2]
        print self.identifingFolder
            
    def FilesInFolder(self):
        files = os.listdir(self.folderPath)
        shots = []
        for i in files:
            if not i.startswith('.~') and i.endswith('.csv'):
                shots.append(i)
        
        # Sort
        timeStep = []
        splice = [8, -4]    #File name lineout<>.csv
        for i in range(len(shots)):
            timeStep.append(int(shots[i][splice[0]:splice[1]]))
        
        timeStep = np.asarray(timeStep)
        sorted_points = sorted(zip(timeStep, shots))
        timeStep = [point[0] for point in sorted_points]
        csv_list = [point[1] for point in sorted_points]
        
        return csv_list, timeStep
        
    def assignSpreadsheetHeader(self):
        if self.mach:
            self.index = ['Ma', 'rho', 'p', 'rho_0', 'Temp', 'p_0', 'U0', 'U_1', 'U_2', 'vtkValidPointMask','arc_length',
             'Points_0','Points_1', 'Points_2']
        else:
            self.index = ['rho', 'p', 'rho_0', 'Temp', 'p_0', 'U0', 'U_1', 'U_2', 'vtkValidPointMask','arc_length',
             'Points_0','Points_1', 'Points_2']
    
    def arrAppend(self, df):
        #Append the correct part of the data to each list
        self.rhoArr.append(df.rho)
        self.pArr.append(df.p)
        self.TArr.append(df.Temp)
        if self.mach:
            self.maArr.append(df.Ma)
        else:
            self.maArr.append(0)
    
    def readInDataToArrs(self):
            #Create lists of the data
        self.rhoArr = []
        self.pArr = []
        self.TArr = []
        self.maArr = []
        
        for csv in self.csv_list:            
            file_interest = self.folderPath + csv
        #    print file_interest
            df = pd.read_csv(file_interest, skiprows=[0], header=None, names=self.index)
            self.arrAppend(df)
            
            #Don't get final data from last file as this may be corrupted if printing crashes.
            if csv is self.csv_list[5]:
                self.arclength = df.arc_length
                self.y_axis = self.arclength * 1e3
        
            
            # First data row is incorrect due to there being no values for the density at time zero
        del self.rhoArr[0]
        del self.maArr[0]
        self.toNpArr()
        
    def toNpArr(self):
        #Convert into numpy arrays
        self.rhoArr = np.asarray(self.rhoArr)
        self.pArr = np.asarray(self.pArr)
        self.TArr = np.asarray(self.TArr)
        self.convert_denity_to_nd()

    def convert_denity_to_nd(self):
        M_A = 28.96e-3
        N_A = 6.022e23
        self.ndArr = (self.rhoArr * N_A) / M_A
        
        
            
    def Evolution_of_ND(self, startingT):
        #A plot showing the lineout evolution with time
        
        fig = plt.figure(figsize=(8, 8)) 
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 2]) 
        
        ax0 = plt.subplot(gs[0])
        plt.title('The evolution of the density profile: ' + self.identifingFolder)
        
        ax1 = plt.subplot(gs[1])
        ax1.set_yscale('log')
        
    
        im = ax0.pcolor(self.ndArr)
        

        #                       [x0, y0, width, height]
        cbar_ax = fig.add_axes([0.805, 0.4, 0.05, 0.45])
        fig.colorbar(im, cax=cbar_ax)
        
    #    fig.colorbar(im, , label=cbarLabel)
    #    cbar_ax = fig.add_axes([1.05, 0.15, 0.05, 0.7])
        
        sumArr = self.ndArr[startingT:]
        FlattenToX = []
        for a in sumArr.T:
            FlattenToX.append(np.average(a))
            
        ax1.plot(self.y_axis, FlattenToX)
        ax0.set_xlim([0, len(FlattenToX)])
        ax1.set_ylabel('Averaged Density for %s to %s' %(startingT, len(self.ndArr)))
        ax0.set_ylabel('Openfoam time dump number')
        ax0.xaxis.set_visible(False)
        ax1.set_xlabel('Lineout Length (mm)')
        xAxis = np.array(bladeline.y_axis)
        ax1.set_xlim([xAxis[0], xAxis[-1]])
        fig.subplots_adjust(right=0.15) 
        plt.tight_layout()
        plt.savefig(self.folderPath + 'Time Evolution of ND ' + self.identifingFolder +'.png', dpi = 250)
        plt.show()
        np.savetxt(self.folderPath + 'Average_lineout_of_shock_' + self.identifingFolder +'.txt' , np.c_[self.y_axis, FlattenToX])
        
        
        indivualLineouts = False
        if indivualLineouts:        
            lineouts = 8
            f, axarr = plt.subplots(nrows = lineouts/2, ncols = 2, figsize=(6, 9))
            
            plotNumbers = np.arange(startingT, len(self.ndArr), (len(self.ndArr) - startingT) / lineouts, dtype = int)
            while len(plotNumbers) > lineouts:
                plotNumbers = plotNumbers[:-1]
            plotNumbers = plotNumbers.reshape(lineouts/2, 2)
            
            
            for row, plot in zip(axarr, plotNumbers):
                for col, p in zip(row, plot):
                    col.plot(self.y_axis, self.ndArr[p])
                    col.title.set_text('TimeStep ' + str(p))
                    col.set_xlabel('Arc length (mm)')
                    np.savetxt(self.folderPath + str(p) + 'lineout.txt', np.c_[self.y_axis, self.ndArr[p]])  
                
            plt.tight_layout()
            plt.suptitle('The evolution of the density profile: ' + self.identifingFolder)
            plt.show()
    
        
mainDir = '/Volumes/CIDU_passport/openFOAM/Line Out Data/PressureScan_5mm_above_Nozzle/'
#mainDir = '/Volumes/CIDU_passport/openFOAM/Line Out Data/1mm_Above_Blade_parallel/'

listSubFolders =  [x[0] for x in os.walk(mainDir)][1:]
for i in range(len(listSubFolders)):
    listSubFolders[i] += '/'
    print listSubFolders[i]
    
#farr = ['h80_0mm/',  'h80_1mm/']
nameSplice = [4, -1]
mach = False
for path in listSubFolders:
    bladeline = LineoutAnalysis(path, mach)
    bladeline.readInDataToArrs()
    bladeline.Evolution_of_ND(5)
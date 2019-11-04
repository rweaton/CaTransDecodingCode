#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 15:38:09 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd
from PeriEventTraceFuncLib import *
from collections import defaultdict
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()

exec(open('./RestoreShelvedWorkspaceScript.py').read())
     
drive, path_and_file = os.path.splitdrive(RestoreFilePath)
path, file = os.path.split(path_and_file)

FigureTitle = file[0:19]
 
# Initialize figure
fig1, axs1 = plt.subplots()
fig1.suptitle(FigureTitle)
 
# Plot performance and performance control plots
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'measure_SE': 'performance_SE',
                'color':'blue'}
 
GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
                     axs1, 'fw_sliding')
 
PlotSpecDict = {'measure': 'performance_median',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'measure_SE': 'performance_SE',
                'color':'lightblue'}
 
GenerateConfIntsPlot(EventsShuffled, EventsShuffled, PlotSpecDict, 
                     axs1, 'fw_sliding')
 
axs1.set_xbound(lower=ParamsDict['BoundaryWindow'][0], 
                upper=ParamsDict['BoundaryWindow'][1])
axs1.set_ybound(lower=0.4, upper=1.)

# Save figure
PerfFigSavePath = (path + os.sep + file[0:19] + '_Perf_' + 'SW_w' + 
                   str(int(1000.*ParamsDict['WindowWidth'])) + 'ms_i' + 
                   str(int(1000.*ParamsDict['StepWidth'])) + 'ms.svg')

fig1.savefig(PerfFigSavePath)

# Plot mutual information and mutual information control plots
fig2, axs2 = plt.subplots()
fig2.suptitle(FigureTitle)
 
PlotSpecDict = {'measure': 'mutual_info',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'measure_SE': 'mutual_info_SE',
                'color':'blue'}
 
GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
                     axs2, 'fw_sliding')
 
PlotSpecDict = {'measure': 'mutual_info_median',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'measure_SE': 'mutual_info_SE',
                'color':'lightblue'}
 
GenerateConfIntsPlot(EventsShuffled, EventsShuffled, PlotSpecDict, 
                     axs2, 'fw_sliding')
 
axs2.set_xbound(lower=ParamsDict['BoundaryWindow'][0], 
                upper=ParamsDict['BoundaryWindow'][1])
axs2.set_ybound(lower=0., upper=1.)

# Save figure
MutInfFigSavePath = (path + os.sep + file[0:19] + '_MutInf_' + 'SW_w' + 
                     str(int(1000.*ParamsDict['WindowWidth'])) + 'ms_i' + 
                     str(int(1000.*ParamsDict['StepWidth'])) + 'ms.svg')

fig2.savefig(MutInfFigSavePath)
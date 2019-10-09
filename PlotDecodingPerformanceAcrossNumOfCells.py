#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 02:12:19 2019

@author: thugwithyoyo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shelve
import os

import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()

exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Determine parent directory and filename from complete path.
drive, path_and_file = os.path.splitdrive(RestoreFilePath)
path, file = os.path.split(path_and_file)

fig1, ax1 = plt.subplots()
fig1.suptitle(FigureTitle)

# Plot decoding on observed activity-outcome correspondence 
PlotSpecDict = {'color':'blue'}
TraceLabel = 'Observed outcomes'
ax1.fill_between(X, PerfFillBand[0,:], PerfFillBand[1,:], 
                        label=TraceLabel, alpha=0.7, color=PlotSpecDict['color'])

#TraceLabel = 'Mean perf.'
ax1.plot(X, PerfMeans, '.-', color=PlotSpecDict['color'])

# Plot decoding performance on activity coupled with shuffled outcomes
TraceLabel = 'Shuffled outcomes'
PlotSpecDict = {'color':'lightblue'}
ax1.fill_between(X, ShuffledPerfFillBand[0,:], ShuffledPerfFillBand[1,:], 
                        label=TraceLabel, alpha=0.7, color=PlotSpecDict['color'])

ax1.plot(X, ShuffledPerfMeans, '.-', color=PlotSpecDict['color'])

# Set plot axes limits and generate labels.
ax1.set_xlabel('number of cells')
ax1.set_ylabel('performance')
ax1.set_ylim([0.4, 1])
ax1.legend(loc='lower right')

# Save figure
fig1.savefig(path + os.sep + file[0:19] + '_PerfVsNumCellsIncluded.svg')
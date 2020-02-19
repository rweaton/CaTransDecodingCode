#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 07:46:53 2019

@author: thugwithyoyo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import CalciumImagingFluorProcessing as CIFP
#import SlidingWindowAnalysisFunc as SWAF
import CalciumTraceDataframeFuncLib as CTDFL
#from PeriEventTraceFuncLib import *
import PeriEventTraceFuncLib as PETFL

from collections import defaultdict
import shelve
import os

import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'

# Have user select workspace to load.  Workspace must be a shelve object '.dat'
# file generated from SlidingWindowAnalysisFunc.py
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()

# Open workspace using Shelve loading script.  Place in try/except block
# to try and bypass occasional loading errors.
#exec(open('./RestoreShelvedWorkspaceScript.py').read())
try:

    exec(open('./RestoreShelvedWorkspaceScript.py').read())
    
except:

    print('Unshelving error.  Will attempt to continue...')

# Determine parent directory and filename from complete path.
drive, path_and_file = os.path.splitdrive(RestoreFilePath)
path, file = os.path.split(path_and_file)

# Define ParamsDict, the dictionary that contains the parameters for
# PLS decoding, Bootstrapping, Shuffle control processes.
# NOTE: ParamsDict is a variable contained in the loaded workspace. Uncomment
# any of the following ParamsDict assignments to make changes to parameter 
# settings for this number-of-cells-dependent analysis.

## Begin settings re-assignment section ## 
#ParamsDict = defaultdict(dict)

# Peripheral target entry events
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
#ParamsDict['AssignedEventVals'] = [-1, 1]

# Set parameters for peri-event extraction
ParamsDict['BoundaryWindow'] = [-1., 1.]
#ParamsDict['StepWidth'] = 0.1
#ParamsDict['WindowWidth'] = 0.4

#ParamsDict['NumLatents'] = 5
#ParamsDict['NumRepetitions'] = 30
#ParamsDict['ConfLevel'] = 0.95
#ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)
## End settings assignment section ##

# Define total number of traces to  be processed per datapoint.
NumTotalTraces = 500

# Load calcium transient data into dataframe.
#CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)

# Load behavioral events into data dictionary
#BehavDict = BehavDictGen(PathToBehavFile)
RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                 'AssignedEventVals' : ParamsDict['AssignedEventVals']}

# Generate an array of dicts with each dict element containing the number of 
# draws to make per datapoint, as well as the number of traces to select per 
# draw.  
# First determine number of cells in cell fluorescence traces dataframe.
# Remember that the first column is a list of timestamps.
NumCells = CellFluorTraces_Frame.shape[1] - 1

# Generate a list array of the number of cell-traces to select per repetition
DrawCounts = np.arange(5, (NumCells + 1), 5)
DrawCounts = np.hstack([np.array([1, 2, 3, 4]), DrawCounts])

# Initialize an array of dictionaries to contain randomly drawn trace samples
(NumSamplingDicts,) = DrawCounts.shape
SamplingDicts = np.empty((NumSamplingDicts,), dtype=dict)

# Populate dictionaries of randomly-drawn traces.  Note that the number of draws
# to perform depends on the number of cells selected.  Here we use the relation:
# NumOfDraws = NumTotalTraces / NumTracesToDraw to determine the number of times 
# to repeat a random draw of a size NumTracesToDraw.  
for i in np.arange(0, NumSamplingDicts):
    
    SamplingDicts[i] = {'NumTracesToDraw': DrawCounts[i],
                        'NumOfDraws' : int(np.round(NumTotalTraces / DrawCounts[i]))}

# Perform decoding on each of the randomly-drawn samples.
for i in np.arange(0, NumSamplingDicts):
        
    # Initialze arrays to contain decoder performance output for each sample.
    PerfDist = np.empty((SamplingDicts[i]['NumOfDraws'],))
    ShuffledPerfDist = np.empty((SamplingDicts[i]['NumOfDraws'],))
    
    # Repeat random sample draws NumOfDraw times.  Each draw selects at random
    # a total of NumTracesToDraw traces.
    for DrawNum in np.arange(0, SamplingDicts[i]['NumOfDraws']):
    
        # Call routine to randomly select traces from cell fluorescence 
        # dataframe.  Output is also a Pandas dataframe.
        SubsetFluorTraces_Frame = CTDFL.TraceSampler(CellFluorTraces_Frame, 
                                    SamplingDicts[i]['NumTracesToDraw'])
        
        # Extract traces from the desired time domain on which decoding will
        # be performed.
        PeriEventExtractorDict = PETFL.PeriEventExtractor_Trace(BehavDict, 
                                    SubsetFluorTraces_Frame, RefEventsDict, 
                                    ParamsDict['BoundaryWindow'])
        
        # Generate a set of indices to test the inclusion portion of the 
        # performance code.
        #InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))
        
        # Determine size of peri-event activity array
        (NumTotalTrials, NumTotalFeatures) = PeriEventExtractorDict['PEA_Array'].shape
        
        # Apply decoding routine to trace snippet.  Store performance output
        # in output array of dicts
        PerfDist[DrawNum] = PETFL.PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    ParamsDict['NumLatents'])['performance']
        
        # Shuffle outcome labels.
        ShuffledIndices = np.arange(0, PeriEventExtractorDict['TargetsVec'].shape[0])
        np.random.shuffle(ShuffledIndices)
        PeriEventExtractorDict['TargetsVec'] = \
            PeriEventExtractorDict['TargetsVec'][ShuffledIndices]
        
        # Rerun the decoder using the list of shuffled outcomes and record
        # performance in the shuffled output array of dicts.  
        ShuffledPerfDist[DrawNum] = PETFL.PLS_DecoderPerformance(
                PeriEventExtractorDict, ParamsDict['NumLatents'])['performance']       

    # Add outcomes arrays for this particular value of draw size to the
    # corresponding array of dicts.
    SamplingDicts[i].update({'PerfDist':PerfDist})
    SamplingDicts[i].update({'ShuffledPerfDist':ShuffledPerfDist})
 
# Extract session id information for figure title.
FigureTitle = file[0:19]

# Initialize arrays to be plotted.  These will contain averages and standard
# error values for both the observed and shuffled datasets.
PerfMeans = np.empty((NumSamplingDicts,))
PerfSEs = np.empty((NumSamplingDicts,))
ShuffledPerfMeans = np.empty((NumSamplingDicts,))
ShuffledPerfSEs = np.empty((NumSamplingDicts,))

# Initialize array to contain corresponding x-axis values for plotting
X = np.empty((NumSamplingDicts,))

# Initialize arrays to contain fill band boundaries for the two trace types
PerfFillBand = np.empty((2, NumSamplingDicts))
ShuffledPerfFillBand = np.empty((2, NumSamplingDicts))

# Iterate through SamplingDicts array and calculate element values of arrays
# to be used as plot input.  Populate corresponding plot arrays that were
# initialized above.
for i in np.arange(0, NumSamplingDicts):

    X[i] = SamplingDicts[i]['NumTracesToDraw']
    PerfMeans[i] = np.mean(SamplingDicts[i]['PerfDist'])
    PerfSEs[i] = (np.std(SamplingDicts[i]['PerfDist'], ddof=1) / 
                  np.sqrt(SamplingDicts[i]['PerfDist'].shape[-1]))
    
    ShuffledPerfMeans[i] = np.mean(SamplingDicts[i]['ShuffledPerfDist'])
    ShuffledPerfSEs[i] = (np.std(SamplingDicts[i]['ShuffledPerfDist'], ddof=1) / 
                  np.sqrt(SamplingDicts[i]['ShuffledPerfDist'].shape[-1]))

# Generate ceiling and floor boundaries of errorband for observed trace
PerfFillBand[0,:] = PerfMeans - PerfSEs
PerfFillBand[1,:] = PerfMeans + PerfSEs

# Generate ceiling and floor boundaries of errorband for shuffled trace
ShuffledPerfFillBand[0,:] = ShuffledPerfMeans - ShuffledPerfSEs
ShuffledPerfFillBand[1,:] = ShuffledPerfMeans + ShuffledPerfSEs

###### Save selected workspace variables ########
# Construct full path to write shelve file
SavePath = path + os.sep + file[0:19] + '_PerfVsNumCellsIncluded'

# Build list of variables to save.
VarsToSave = ['SamplingDicts', 'PerfMeans', 'ShuffledPerfMeans', 
              'PerfSEs', 'ShuffledPerfSEs', 'X', 'FigureTitle', 
              'PerfFillBand', 'ShuffledPerfFillBand', 'ParamsDict', 
              'NumTotalTraces', 'NumCells', 'RefEventsDict']

# Initialize and open self object
my_shelf = shelve.open(SavePath)

# Iterate through the list of keys in the save list, writing each to shelf
# object
for key in VarsToSave:
    
    # Surround write command in try/except clause to bypass possible TypeError
    # issues.
    try:
        
        #my_shelf[key] = globals()[key]
        my_shelf[key] = locals()[key]
        
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
 
# Close shelf object after variables have been written to file
my_shelf.close()
###### End workspace save ###########

#### Begin line plot generation #####
fig1, ax1 = plt.subplots()
fig1.suptitle(FigureTitle)

## Plot decoding of observed activity-outcome correspondence ##

# Set plot color of observed trace
PlotSpecDict = {'color':'orange'}

# Label observed trace for legend
TraceLabel = 'Observed outcomes'

# 
ax1.fill_between(X, PerfFillBand[0,:], PerfFillBand[1,:], 
                        label=TraceLabel, alpha=0.7, color=PlotSpecDict['color'])

#TraceLabel = 'Mean perf.'

ax1.plot(X, PerfMeans, '.-', color=PlotSpecDict['color'])

## Plot decoding performance on activity coupled with shuffled outcomes ##

# Set plot color of shuffled trace
PlotSpecDict = {'color':'gray'}

# Label shuffled trace for legend
TraceLabel = 'Shuffled outcomes'

# Plot errorband for shuffled trace
ax1.fill_between(X, ShuffledPerfFillBand[0,:], ShuffledPerfFillBand[1,:], 
                        label=TraceLabel, alpha=0.7, color=PlotSpecDict['color'])

# Plot shuffled trace
ax1.plot(X, ShuffledPerfMeans, '.-', color=PlotSpecDict['color'])

# Set plot axes limits and generate labels.
ax1.set_xlabel('number of cells')
ax1.set_ylabel('performance')
ax1.set_ylim([0.4, 1])
ax1.legend(loc='lower right')
####### End line plot generation ###########

####### Save figure ###########
fig1.savefig(path + os.sep + file[0:19] + '_PerfVsNumCellsIncluded.svg')
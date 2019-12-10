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
from PeriEventTraceFuncLib import *
from collections import defaultdict
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

# Define paths to JSON files to load data.
#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-19/2018-12-19-11-24-58_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-19/2018-12-19-11-24-58_new_unique_C.json'

# Define path to save analysis results.
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-19-11-24-58_new_unique_400ms_SamFiltered'

# Define ParamsDict, the dictionary that contains the parameters for
# PLS decoding, Bootstrapping, Shuffle control processes.
# Peripheral target entry events

#ParamsDict = defaultdict(dict)
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
#ParamsDict['AssignedEventVals'] = [-1, 1]

# Set parameters for peri-event extraction
ParamsDict['BoundaryWindow'] = [0., 1.]
#ParamsDict['StepWidth'] = 0.1
#ParamsDict['WindowWidth'] = 0.4

#ParamsDict['NumLatents'] = 5
#ParamsDict['NumRepetitions'] = 30
#ParamsDict['ConfLevel'] = 0.95
#ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)

# Define total number of traces to  be processed per datapoint.
NumTotalTraces = 500

# Load calcium transient data into dataframe.
#CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)

# Load behavioral events into data dictionary
#BehavDict = BehavDictGen(PathToBehavFile)

RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                 'AssignedEventVals' : ParamsDict['AssignedEventVals']}

# Generate an array of dicts with each dict containing the number of draws to
# make per datapoint as well as the number of samples to take per draw
NumCells = CellFluorTraces_Frame.shape[1] - 1

DrawCounts = np.arange(5, (NumCells + 1), 5)
DrawCounts = np.hstack([np.array([1, 2, 3, 4]), DrawCounts])

(NumSamplingDicts,) = DrawCounts.shape
SamplingDicts = np.empty((NumSamplingDicts,), dtype=dict)


for i in np.arange(0, NumSamplingDicts):
    
    SamplingDicts[i] = {'NumTracesToDraw': DrawCounts[i],
                        'NumOfDraws' : int(np.round(NumTotalTraces / DrawCounts[i]))}


for i in np.arange(0, NumSamplingDicts):
        
    PerfDist = np.empty((SamplingDicts[i]['NumOfDraws'],))
    ShuffledPerfDist = np.empty((SamplingDicts[i]['NumOfDraws'],))
    
    for DrawCount in np.arange(0, SamplingDicts[i]['NumOfDraws']):
    
        SubsetFluorTraces_Frame = CTDFL.TraceSampler(CellFluorTraces_Frame, 
                                    SamplingDicts[i]['NumTracesToDraw'])

        PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                        SubsetFluorTraces_Frame, RefEventsDict, 
                                        ParamsDict['BoundaryWindow'])
        
        # Generate a set of indices to test the inclusion portion of the performance code.
        (NumTotalTrials, NumTotalFeatures) = PeriEventExtractorDict['PEA_Array'].shape
        #InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))
    
        PerfDist[DrawCount] = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    ParamsDict['NumLatents'])['performance']
        
        ShuffledIndices = np.arange(0, PeriEventExtractorDict['TargetsVec'].shape[0])
        np.random.shuffle(ShuffledIndices)
        PeriEventExtractorDict['TargetsVec'] = \
            PeriEventExtractorDict['TargetsVec'][ShuffledIndices]
        
        ShuffledPerfDist[DrawCount] = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    ParamsDict['NumLatents'])['performance']       

    SamplingDicts[i].update({'PerfDist':PerfDist})
    SamplingDicts[i].update({'ShuffledPerfDist':ShuffledPerfDist})
 
# Begin plotting
FigureTitle = file[0:19]
    
PerfMeans = np.empty((NumSamplingDicts,))
PerfSEs = np.empty((NumSamplingDicts,))
ShuffledPerfMeans = np.empty((NumSamplingDicts,))
ShuffledPerfSEs = np.empty((NumSamplingDicts,))
X = np.empty((NumSamplingDicts,))

PerfFillBand = np.empty((2, NumSamplingDicts))
ShuffledPerfFillBand = np.empty((2, NumSamplingDicts))

for i in np.arange(0, NumSamplingDicts):

    X[i] = SamplingDicts[i]['NumTracesToDraw']
    PerfMeans[i] = np.mean(SamplingDicts[i]['PerfDist'])
    PerfSEs[i] = (np.std(SamplingDicts[i]['PerfDist'], ddof=1) / 
                  np.sqrt(SamplingDicts[i]['PerfDist'].shape[-1]))
    
    ShuffledPerfMeans[i] = np.mean(SamplingDicts[i]['ShuffledPerfDist'])
    ShuffledPerfSEs[i] = (np.std(SamplingDicts[i]['ShuffledPerfDist'], ddof=1) / 
                  np.sqrt(SamplingDicts[i]['ShuffledPerfDist'].shape[-1]))

PerfFillBand[0,:] = PerfMeans - PerfSEs
PerfFillBand[1,:] = PerfMeans + PerfSEs

ShuffledPerfFillBand[0,:] = ShuffledPerfMeans - ShuffledPerfSEs
ShuffledPerfFillBand[1,:] = ShuffledPerfMeans + ShuffledPerfSEs

# Save selected workspace variables
SavePath = path + os.sep + file[0:19] + '_PerfVsNumCellsIncluded'
VarsToSave = ['SamplingDicts', 'PerfMeans', 'ShuffledPerfMeans', 
              'PerfSEs', 'ShuffledPerfSEs', 'X', 'FigureTitle', 
              'PerfFillBand', 'ShuffledPerfFillBand', 'ParamsDict', 
              'NumTotalTraces', 'NumCells', 'RefEventsDict']

my_shelf = shelve.open(SavePath)

for key in VarsToSave:
    
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
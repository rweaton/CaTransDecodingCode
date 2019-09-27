#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 23:40:34 2019

@author: thugwithyoyo
"""

import matplotlib.pyplot as plt
from PeriEventTraceFuncLib import *
import numpy as np
from CaTraceNormalizer import *
from collections import defaultdict

RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56_unique_800ms_SlidingWindow-filtered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56/SW_400ms_100ms/2018-12-10-11-37-56_unique_400ms_SlidingWindow.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-14-11-01-41/JointHemisphereDecoding/SingleArmDirectionTask/2018-12-14-11-01-41_new_unique_SW_LHem.dat'
exec(open('./RestoreShelvedWorkspaceScript.py').read())

ParamsDict['SearchDomain'] = [-1., 1.]
ColorMap = 'seismic'

SortingRefEvent = RefEventsList[1]

RefEventsList = RefEventsDict['RefEventsList']
nEvents = len(RefEventsList)

# Acquire number of columns in fluorescence trace dataframe
nCols = CellFluorTraces_Frame.shape[1]

# zScore all calcium traces prior to peri-event extraction
FullTraces = pd.DataFrame(CellFluorTraces_Frame._slice(slice(1, nCols), 1)).values.transpose()

zScoredTraces = zScoreTraces(FullTraces, ParamsDict).transpose()

# Replace raw fluorescence traces with z-scores in the trace dataframe
IndicesOfColsToChange = np.arange(1, nCols)
ColNames = list(CellFluorTraces_Frame.columns)

for i in IndicesOfColsToChange:
    CellFluorTraces_Frame[ColNames[i]] = zScoredTraces[:, i - 1]
    
#CellFluorTraces_Frame._slice(slice(1, nCols), 1).values = zScoredTraces.transpose() 

# Recompile PeriEvent activity dict to include Full Domain in boundary window.
PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                                  CellFluorTraces_Frame, 
                                                  RefEventsDict, 
                                                  ParamsDict['BoundaryWindow'])

NumTrials, NumCellsByNumSamples = PeriEventExtractorDict['PEA_Array'].shape

(_, NumColumns) = CellFluorTraces_Frame.shape

NumCells = NumColumns - 1

NumSamples = int(NumCellsByNumSamples/NumCells)

RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], ParamsDict['BoundaryWindow'][1], num=NumSamples, endpoint=False)

NumEvents = len(RefEventsList)

AveragedTracesMatrices = np.empty((NumCells, NumSamples, NumEvents))
PeakNormalizedTraces = np.empty((NumCells, NumSamples, NumEvents))

fig, axes  = plt.subplots(1, NumEvents)

# generate 2 2d grids for the x & y bounds
xv, yv = np.meshgrid(RelTimeVec, np.arange(1, NumCells + 1, 1))
i = 0

SortProcessingDict = defaultdict(dict)

for RefEvent in RefEventsList:
    
    TracesFilt = PeriEventExtractorDict['TrialIndicesByEventDict'][RefEvent]
    
    AveragedTraces = np.mean(PeriEventExtractorDict['PEA_Array'][TracesFilt,:],
                                      axis=0)
    

    
    AveragedTracesMatrices[:,:,i] = np.reshape(AveragedTraces, (NumCells, NumSamples))
    
    
    SortProcessingDict[RefEvent] = CaTraceSorter(AveragedTracesMatrices[:,:,i],
                                                 ParamsDict, 'PeakSort')
    
    PeakNormalizedTraces[:, :, i] = PeakNormalizer(AveragedTracesMatrices[:,:,i])

#    SortProcessingDict[RefEvent] = CaTraceSorter(PeakNormalizedTraces[:,:,i],
#                                                 ParamsDict, 'ActivitySort')
    
    i += 1
    
#vmin = np.min(np.min(np.min(AveragedTracesMatrices, axis=0),axis=0),axis=0)
# Determine maximum value to set color map.  Take floor to the next 0.5 level
vmax = np.max(np.max(np.max(AveragedTracesMatrices, axis=0),axis=0),axis=0)
vmax = np.floor(2.0*vmax)/2.0
vmin = -vmax

#vmin = np.min(np.min(np.min(PeakNormalizedTraces, axis=0),axis=0),axis=0)
#vmax = np.max(np.max(np.max(PeakNormalizedTraces, axis=0),axis=0),axis=0)

i = 0
for RefEvent in RefEventsList:
    
    
    im = axes[i].pcolor(xv, yv, AveragedTracesMatrices[SortProcessingDict[SortingRefEvent]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)
#    im = axes[i].pcolormesh(xv, yv, PeakNormalizedTraces[SortProcessingDict[RefEventsList[0]]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap='coolwarm')

    #fig.colorbar(im, ax=axes[i])
    axes[i].set_xlabel('time relative to reach (sec.)')
    
    if i == 0:
        axes[0].set_ylabel('Cell identity')

#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']    
    if (RefEvent == 'M6T0_Entry_ts'):
        SubPlotTitle = 'RH_Zone1'
    elif(RefEvent == 'M6T1_Entry_ts'):
        SubPlotTitle = 'RH_Zone2'
    elif(RefEvent == 'M7T0_Entry_ts'):
        SubPlotTitle = 'LH_Zone1'
    elif(RefEvent == 'M7T1_Entry_ts'):
        SubPlotTitle == 'LH_Zone2'
        
    axes[i].title.set_text(SubPlotTitle)
    
    i += 1

fig.colorbar(im, ax=list(axes), orientation='horizontal')

# Define figure name
drive, path_and_file = os.path.splitdrive(RestoreFilePath)
path, file = os.path.split(path_and_file)
FigureTitle = file[0:19] + ' trial-averaged z-scores of $\Delta$F/F Ca response'
fig.suptitle(FigureTitle)

#  Generate difference plot. 
#fig2, axs2 = plt.subplots(1,1)
#
#DiffMat = np.diff(AveragedTracesMatrices, axis=-1)[:,:,0]
##DiffMat = np.diff(PeakNormalizedTraces, axis=-1)[:,:,0]
#
#im2 = axs2.pcolormesh(xv, yv, DiffMat[SortProcessingDict[SortingRefEvent]['SortIndices'],:], cmap=ColorMap)
#fig2.colorbar(im2, ax=axs2)
#axs2.set_ylabel('Cell identity')
#axs2.set_xlabel('time relative to reach (sec.)')

# Generate a tuning plot
AbsSortProcessingDict = defaultdict(dict)

NormVals = np.zeros((NumCells,))

i = 0
for RefEvent in RefEventsList:
    
    AbsSortProcessingDict[RefEvent] = CaTraceSorter(AveragedTracesMatrices[:,:,i],
                                                    ParamsDict, 'AbsPeakSort')
    
    NormVals = NormVals + np.abs(AbsSortProcessingDict[RefEvent]['MaxVals'])
    
    i += 1

NormVals = NormVals / len(RefEventsList)

TuningMatrix = np.divide(np.diff(AveragedTracesMatrices, axis=-1)[:,:,0],
                         np.dot(np.array([NormVals]).transpose(), 
                         np.array([np.ones(NumSamples,)])))

vmax = np.max(np.max(TuningMatrix, axis=0),axis=0)
vmax = np.floor(2.0*vmax)/2.0
vmin = -vmax

fig2, axs2 = plt.subplots(1,1)
im2 = axs2.pcolormesh(xv, yv, TuningMatrix[SortProcessingDict[SortingRefEvent]['SortIndices'],:], cmap='PiYG', vmin=vmin, vmax=vmax)
fig2.colorbar(im2, ax=axs2)
axs2.set_ylabel('Cell identity')
axs2.set_xlabel('time relative to reach (sec.)')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:27:34 2019

@author: thugwithyoyo
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from PeriEventTraceFuncLib import *
import numpy as np
from CaTraceNormalizer import *
from collections import defaultdict

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import messagebox as mbox 


def GetProportionsOfTunedCells():
    
    root = tk.Tk()
    RestoreFilePath = askopenfilename()
    root.withdraw()
    
    # Determine parent directory and filename from complete path.
    Drive, Path_and_file = os.path.splitdrive(RestoreFilePath)
    Path, File = os.path.split(Path_and_file)
    
    # Attempt to load decoding output workspace into function scope.
    try:
    
        exec(open('./RestoreShelvedWorkspaceScript.py').read())

    except:
        
        print('Error when restoring shelved workspace.  Will attempt to continue')
    
    ### Specify a dictionary of parameters for tuning calcluation and plotting.
    ParamsDict = defaultdict(dict)
    #ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
    ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']
    
    # Scalar values assigned to event types listed above.
    ParamsDict['AssignedEventVals'] = [-1, 1]
    
    # Set parameters for peri-event extraction
    ParamsDict['BoundaryWindow'] = [-1., 2.]
    ParamsDict['StepWidth'] = 0.1
    ParamsDict['WindowWidth'] = 0.4
    
    ParamsDict['NumLatents'] = 5
    ParamsDict['NumRepetitions'] = 30
    ParamsDict['ConfLevel'] = 0.95
    ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)
    
    # Time domain to consider for calculation of tuning index
    ParamsDict['SearchDomain'] = [-1., 1.]
    
    # Specify the method to be used to calculate tuning method. 
    #ParamsDict['TuningType'] = 'PlainDifference'
    #ParamsDict['TuningType'] = 'WeightedDifference'
    ParamsDict['TuningType'] = 'DiffOfAvgsOverSumOfMagOfAvgs'
    
    # For the WeightedDifference method, specify the standard deviation of the 
    # Gaussian weighting function.
    ParamsDict['StdDev'] = 0.75
    
    # Specify the number of bins to use to generate the histogram.
    ParamsDict['NumHistoBins'] = 20
    
    # Specify the list of quantiles to plot as vertical lines on histogram
    ParamsDict['QuantilesList'] = [0.25, 0.75]
    
    # The magnitude of the tuning index value that separates tuned cells from 
    # non-tuned cells.
    ParamsDict['TuningCutoffLevel'] = 0.75
    #ParamsDict['TuningCutoffLevel'] = 0.2   

    RefEventsList = RefEventsDict['RefEventsList']
    SortingRefEvent = RefEventsList[1]
    
    nEvents = len(RefEventsList)
    
    # Load in CentroidFrame just for the 2018-12-11-... analysis. Shelve is wonkey.
    # !!!! Be sure to comment this out after!!!
    #CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)
    
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
    
    # Count the number of trials and determine length of each row in the formatted
    # Peri-event activity array.
    NumTrials, NumCellsByNumSamples = PeriEventExtractorDict['PEA_Array'].shape
    
    # Count the total number of columns in the fluorescence trace dataframe
    (_, NumColumns) = CellFluorTraces_Frame.shape
    
    # 
    NumCells = NumColumns - 1
    
    NumSamples = int(NumCellsByNumSamples/NumCells)
    
    RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
                             ParamsDict['BoundaryWindow'][1], 
                             num=NumSamples, endpoint=False)
    
    NumEvents = len(RefEventsList)
    
    AveragedTracesMatrices = np.empty((NumCells, NumSamples, NumEvents))
        
    SortProcessingDict = defaultdict(dict)
    
    i = 0
    for RefEvent in RefEventsList:
        
        TracesFilt = PeriEventExtractorDict['TrialIndicesByEventDict'][RefEvent]
        
        AveragedTraces = np.mean(PeriEventExtractorDict['PEA_Array'][TracesFilt,:],
                                          axis=0)
            
        AveragedTracesMatrices[:,:,i] = np.reshape(AveragedTraces, (NumCells, NumSamples))
        
        i += 1
   
    TuningScalars = GetTuning(AveragedTracesMatrices[:,:,::-1], 
                          ParamsDict['TuningType'], 
                          ParamsDict=ParamsDict) 
               
    SortProcessingDict['AvgVals'] = np.array([TuningScalars]).transpose()
    SortProcessingDict['SortIndices'] = np.argsort(SortProcessingDict['AvgVals'], axis=0).transpose()[0]
        
    ###################################
    Proportions = np.empty((3,))
    
    (TotalNumIndices,_) = SortProcessingDict['AvgVals'].shape
    
    # Calulate the proportion of Zone2-tuned cells
    Proportions[0] = np.sum(SortProcessingDict['AvgVals'] <=
               -1.*ParamsDict['TuningCutoffLevel'])/TotalNumIndices
    
    # Calculate the proportion of Zone1-tuned cells
    Proportions[2] = np.sum(SortProcessingDict['AvgVals'] >=  
               ParamsDict['TuningCutoffLevel'])/TotalNumIndices
    
    # Calculate the remaining proportion  of non-selective cells.
    Proportions[1] = 1. - (Proportions[0] + Proportions[2])
    
    SessionName = File[0:19]
    
    return Proportions, SessionName

# Debug run
#Proportions, SessionName = GetProportionsOfTunedCells(ParamsDict)

# Have user build array of proportions across sessions by repeatedly selecting
# sessions to include and appending them to group arrays.

ProportionsArray = np.array([], dtype=float)
SessionNamesArray = np.array([], dtype=object)

AddAnotherSession = True

while AddAnotherSession == True:
    
    ProportionsVec, SessionName = GetProportionsOfTunedCells()
    
    if ProportionsArray.shape == (0,):
        
        ProportionsArray = ProportionsVec
        SessionNamesArray = SessionName
        
    else:   
        
        ProportionsArray = np.vstack((ProportionsArray, ProportionsVec))
        SessionNamesArray = np.hstack((SessionNamesArray, SessionName))

    root = tk.Tk()
     
    AddAnotherSession = mbox.askyesno('Add session user query', 'Add another session?')
    
    root.withdraw()

SortMap = np.argsort(SessionNamesArray)
SessionNamesArray = SessionNamesArray[SortMap]
ProportionsArray = ProportionsArray[SortMap,:]
    
fig, axs = plt.subplots(1,1) 

PlotType = 'errorbars'
(NumDatapoints, NumCategories) = ProportionsArray.shape
xVals = np.outer(np.ones_like(ProportionsArray[:,0]),np.array([[1, 2, 3]]))
Labels = np.array(['zone2 sel.', 'non-tuned', 'zone1 sel.'])
ScatterDatapointsLabel='session prop.'
ColumnOrder = np.array([1, 2, 0])
SymbolSet = np.array(['o', 'v', '^', 's', 'p', 'P', '*', 'H', 'X', 'D', 'd'])

# Scatter plot visual parameter settings
MarkerSize = 80

# Barplot visual parameter settings
ColorsSet = np.array(['green', 'lightgray', 'magenta'])

# Boxplot visual parameters settings
boxprops = dict(linestyle='-', linewidth=3, color='blue')

capandwhiskerprops = dict(linestyle='-', linewidth=1, color='black')

flierprops = dict(marker='o', markerfacecolor='black', markersize=6,
                  linestyle='none', alpha=0.5)

medianprops = dict(linestyle='-', linewidth=3, color='red')

meanpointprops = dict(marker='D', markeredgecolor='black',
                      markerfacecolor='firebrick')
meanlineprops = dict(linestyle='--', linewidth=2.5, color='purple')
    
if PlotType == 'boxplot':
    
    for i in np.arange(0, NumDatapoints):
        
        axs.scatter(xVals[i, :], ProportionsArray[i, ColumnOrder], 
                    marker=SymbolSet[i], color='gray', s=MarkerSize, alpha=0.5, 
                    label=SessionNamesArray[i])
    
        bp = axs.boxplot(ProportionsArray[:, ColumnOrder], 
                         labels=Labels[ColumnOrder], notch=False, showfliers=False, 
                         bootstrap=1000, boxprops=boxprops, flierprops=flierprops, 
                         medianprops=medianprops, whiskerprops=capandwhiskerprops, 
                         capprops=capandwhiskerprops)
    
elif PlotType == 'errorbars': 

    PropMeans = np.mean(ProportionsArray, axis=0)
    StdDevs = np.std(ProportionsArray, axis=0)
    StdErrors = StdDevs/np.sqrt(NumDatapoints)
    
    barlist = axs.bar(xVals[0, :],  PropMeans[ColumnOrder], alpha=0.5)
    
    for j in np.arange(0, NumCategories):
        
        barlist[j].set_color(ColorsSet[ColumnOrder][j])
    
    for i in np.arange(0, NumDatapoints):
        
            axs.scatter(xVals[i, :], ProportionsArray[i, ColumnOrder], 
                    marker=SymbolSet[i], color='gray', s=MarkerSize, alpha=0.5, 
                    label=SessionNamesArray[i])
    
    axs.errorbar(xVals[0, :],  PropMeans[ColumnOrder], yerr=StdDevs, fmt='.',
                 color='black', ecolor='black')
    

# Set x tick labels
axs.set_xticks(xVals[0, :])
axs.set_xticklabels(Labels[ColumnOrder])

# Remove bounding box
axs.spines['right'].set_visible(False)
axs.spines['top'].set_visible(False)

axs.set_ylabel('proportion of cells')
axs.set_ylim([0., 1.])

# Add legend
axs.legend(prop={'size': 8})

# Add figure title
fig.suptitle('Session distributions of cell direction preferences')

################ Begin old version ################

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-11-21/2018-11-21-10-49-56_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-11-21/2018-11-21-10-49-56_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-11-21-10-49-56/2018-11-21-10-49-56_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-11-26/2018-11-26-11-45-46_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-11-26/2018-11-26-11-45-46_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-11-26-11-45-46/2018-11-26-11-45-46_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-04/2018-12-04-10-31-21_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-04/2018-12-04-10-31-21_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-04-10-31-21/2018-12-04-10-31-21_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10/2018-12-10-11-37-56_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10/2018-12-10-11-37-56_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56/2018-12-10-11-37-56_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-11/2018-12-11-10-53-04_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-11/2018-12-11-10-53-04_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-11-10-53-04_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-17/2018-12-17-11-38-42_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-17/2018-12-17-11-38-42_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-17-11-38-42/2018-12-17-11-38-42_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-18/2018-12-18-11-20-21_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-18/2018-12-18-11-20-21_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-18-11-20-21/2018-12-18-11-20-21_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-19/2018-12-19-11-24-58_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-19/2018-12-19-11-24-58_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-19-11-24-58/2018-12-19-11-24-58_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-31-41_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-31-41_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-07-10-31-41/2019-01-07-10-31-41_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-45-52_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-45-52_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-07-10-45-52/2019-01-07-10-45-52_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-24-11-50-23/2019-01-24-11-50-23_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-24-11-36-02/2019-01-24-11-36-02_new_unique_400ms_SamFiltered'

#RecordingLabels = ['2018-11-21-10-49-56',
#                   '2018-11-26-11-45-46',
#                   '2018-12-04-10-31-21',
#                   '2018-12-10-11-37-56',
#                   '2018-12-11-10-53-04',
#                   '2018-12-17-11-38-42',
#                   '2018-12-18-11-20-21',
#                   '2018-12-19-11-24-58',
#                   '2019-01-07-10-31-41',
#                   '2019-01-07-10-45-52',
#                   '2019-01-24-11-36-02',
#                   '2019-01-24-11-50-23',
#                   ]
#NonTunedProps = np.array([45.7, 41.1, 55.7, 49.2, 61.8, 60.2, 58.7, 55.3, 57.3, 50.6, 43.4, 34.9])
#Zone1TunedProps = np.array([28.6, 35.6, 25.2, 20.8, 26.4, 17.9, 15.4, 26.8, 24.7, 24.7, 31.1, 32.1])
#Zone2TunedProps = np.array([25.7, 23.3, 19.1, 30.0, 11.8, 22.0, 26.0, 17.9, 18.0, 24.7, 25.5, 33.0])
#
#fig, axs = plt.subplots(1,1)
#plt.plot(RecordingLabels, NonTunedProps, '.-', label='non-tuned')
#plt.plot(RecordingLabels, Zone1TunedProps, '.-', label='z2 pref.')
#plt.plot(RecordingLabels, Zone2TunedProps, '.-', label='z1 pref.')
#
#axs.set_xlabel('recording ID')
#axs.set_ylabel('percentage of cells')
#plt.xticks(rotation=90)


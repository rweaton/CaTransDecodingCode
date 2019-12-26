#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:33:08 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd
from PeriEventTraceFuncLib import *
from collections import defaultdict

import os

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

# Paths to data in JSON formatted files.  
# NOTE: These have been commented out so that the SlidingWindowAnalysisFunc 
# function can be called without these global level variables being 
# intitialized.

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-31-41_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-31-41_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-07-10-31-41/2019-01-07-10-31-41_new_unique_400ms_SamFiltered'


#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-45-52_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-07/2019-01-07-10-45-52_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-07-10-45-52/2019-01-07-10-45-52_new_unique_400ms_SamFiltered'


################### Analyses run ###############
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
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-11-10-53-04/2018-12-11-10-53-04_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14/2018-12-14-11-01-41_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14/2018-12-14-11-01-41_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-14-11-01-41/2018-12-14-11-01-41_new_unique_400ms_SamFiltered'


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

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-24-11-36-02/2019-01-24-11-36-02_new_unique_400ms_SamFiltered'

#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_C.json'
#SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2019-01-24-11-50-23/2019-01-24-11-50-23_new_unique_400ms_SamFiltered'


# Define ParamsDict, the dictionary that contains the parameters for
# PLS decoding, Bootstrapping, Shuffle control processes.
# Peripheral target entry events

ParamsDict = defaultdict(dict)
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']
ParamsDict['RefEventsList'] = ['M7T0_Entry_ts', 'M7T1_Entry_ts']

# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set parameters for peri-event extraction
ParamsDict['BoundaryWindow'] = [-2., 2.]
ParamsDict['StepWidth'] = 0.1
ParamsDict['WindowWidth'] = 0.4

ParamsDict['NumLatents'] = 5
ParamsDict['NumRepetitions'] = 30
ParamsDict['ConfLevel'] = 0.95
ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)

def SlidingWindowAnalysisFunc(BehavDict, CellFluorTraces_Frame, SavePath, ParamsDict):

    # Peripheral target entry events
    #Params.RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']
    
    # Scalar values assigned to event types listed above.
    #ParamsDict.AssignedEventVals = [-1, 1]
    
    # Pack into a dict the two lists above that specify info for reference 
    # events
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
        
    # Set parameters for peri-event extraction
    #ParamsDict.BoundaryWindow = [-1., 2.]
    #ParamsDict.StepWidth = 0.1
    #ParamsDict.WindowWidth = 0.4
    
    ArrayOfSlidWindows = SlidingWindowGen(ParamsDict['BoundaryWindow'], 
                                          ParamsDict['StepWidth'], 
                                          ParamsDict['WindowWidth'])
    
    # Set parameters for PLS
    #ParamsDict.NumLatents = 5
    
    # Set parameters for Monte Carlo estimation of  confidence intervals
    #ParamsDict.NumRepetitions = 30
    #ParamsDict.ConfLevel = 0.95
    
    # Specified anti-tolerance window, relative to target entry, for detecting and
    # removing repeat entries that followed shortly after the initial entry.
    #ParamsDict.RelativeTolWindow = (0.0001, 2.5)
    
    # Generate the unfiltered behavior dictionary.
    if (type(BehavDict) == str):
        
        PathToBehavFile = BehavDict
        BehavDict = BehavDictGen(PathToBehavFile)
        
    elif (type(BehavDict) == dict):
        pass
    
#    # Detect rapid repeats within each event list.
#    EventFilters = RemoveRepeatTargetEntries(BehavDict, 
#                                             RefEventsDict['RefEventsList'], 
#                                             ParamsDict['RelativeTolWindow'])
#    
#    # Remove repeat events
#    for ef in EventFilters:
#        
#        BehavDict[ef] = BehavDict[ef][EventFilters[ef]]
#        
#    # Determine lowest number of occurances among event types.
#    MinNumOccurences = BehavDict[ParamsDict['RefEventsList'][0]].shape[0]
#    MinEventType = ParamsDict['RefEventsList'][0]
#
#    for e in ParamsDict['RefEventsList']:
#        
#        NumOccurences = BehavDict[e].shape[0]
#        
#        if NumOccurences < MinNumOccurences:
#            
#            MinNumOccurences = NumOccurences
#            MinEventType = e
#            
#    # Randomly draw, without replacement, the equivalent number of event occurances
#    # from longer lists belonging to other event types.
#    Filt = (np.array(ParamsDict['RefEventsList']) != MinEventType)
#    EventListsToSample = np.array(ParamsDict['RefEventsList'])[Filt]
#
#    for e in EventListsToSample:
#        
#        BehavDict[e] = pd.Series(np.random.choice(BehavDict[e], 
#                  MinNumOccurences, replace=False))
    
    # If a path is given in in place of a datafrom for the second arg, then 
    # import the dataframe from filepath string.  
    if type(CellFluorTraces_Frame) == str:
        
        # Generate the data frame of calcium transients.
        PathToFluorFile = CellFluorTraces_Frame
        CellFluorTraces_Frame = CellFluorTraces_FrameGen(PathToFluorFile)
            
        # Assemble dataframe of centroids.
        CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)
        
    elif type(CellFluorTraces_Frame) == pd.DataFrame:
        pass
    
    # Grow window forwards from floor
    (NumDomains, _) = ArrayOfSlidWindows.shape
    
    # Initialize an empty array to contain output dictionaries from the 
    # decoder cross-validation perfomance and monte carlo bootstrap routines
    Performance = np.empty((NumDomains,), dtype=dict)
    ConfInts = np.empty((NumDomains,), dtype=dict)
    EventsShuffled = np.empty((NumDomains,), dtype=dict)
    
    for i in np.arange(0, NumDomains):
     
        PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                        CellFluorTraces_Frame, RefEventsDict, 
                                        ArrayOfSlidWindows[i])
        
        # Generate a set of indices to test the inclusion portion of the performance code.
        PEA_Array = PeriEventExtractorDict['PEA_Array']
        (NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
        #InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))
    
        Performance[i] = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                                ParamsDict['NumLatents'])
        Performance[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
    
        ConfInts[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                     ParamsDict['NumLatents'], 
                                     ParamsDict['NumRepetitions'], 
                                     ParamsDict['ConfLevel'],
                                     Performance[i])
        
        ConfInts[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
        
        EventsShuffled[i] = PLS_Shuffle(PeriEventExtractorDict, 
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        
#        EventsShuffled[i] = PLS_Shuffle2(PeriEventExtractorDict, 
#                                        ParamsDict['NumLatents'], 
#                                        ParamsDict['NumRepetitions'], 
#                                        ParamsDict['ConfLevel'])
                        
        EventsShuffled[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
        
    #####################################
    ########   Start shelving   #########
    #####################################
    
    # Assemble save path
    # get root directory of save path from path to calcium data
    # SavePath is an arguement above.  The following script requires it 
    exec(open('./ShelveWorkspaceScript.py').read())
    
###############################################################################
#########                   START EXECUTION SCRIPT                    #########
###############################################################################
# Acquire path of workspace to load.
# Set defualt parent directories
DataRootDir = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData'
SaveRootDir = '/home/thugwithyoyo/CaTransDecoding/Output'

# Start file dialogs 
root = tk.Tk()

# Prompt user to navigate to, and select, the session behavior file
PathToBehavFile = askopenfilename(title='Select session behavior file',
                                  filetypes=[("json files","*.json"), 
                                             ("csv files", "*.csv")],
                                  initialdir=DataRootDir)

# Prompt user to navigate to, and select, the session imaging file
PathToFluorFile = askopenfilename(title='Select session fluorescence record file',
                                  filetypes=[("json files","*.json"),
                                             ("csv files", "*.csv")],
                                  initialdir=DataRootDir)


# Determine parent directory and filename from complete path.
Drive, Path_and_file = os.path.splitdrive(PathToBehavFile)
Path, File = os.path.split(Path_and_file)

# Extract session id from filename
SessionID = File[0:19]

# Construct a  workspace filename to use as save default
SaveFilename =  SessionID + '_new_unique_400ms_SamFiltered'
DefaultSavePath = SaveRootDir + os.path.sep + SessionID 

# Prompt user to select default
SavePath = asksaveasfilename(title='Set workspace save path',
                             initialdir=DefaultSavePath,
                             initialfile=SaveFilename)

# End file dialogs
root.withdraw()

# Run sliding window analysis using above-defined subroutine.
SlidingWindowAnalysisFunc(PathToBehavFile, PathToFluorFile, SavePath, ParamsDict)
#SlidingWindowAnalysisFunc(PathToBehavFile, FluorDataframe_Combined, SavePath, ParamsDict)
#SlidingWindowAnalysisFunc(BehavDict_Combined, CaImag_df_Combined, SavePath, ParamsDict)

# Figure generation code below moved to GenerateDecodePlotsScript.py.  Use the
# script to generate plots from shelved workspaces generated from above-called function.
# 
#RestoreFilePath = SavePath +'.dat'
#exec(open('./RestoreShelvedWorkspaceScript.py').read())
#     
##### Plot outcome measures #####
# 
## Specify outcome measures to  be plotted
##PerformancePlotSpecDict = {'measure': 'performance',
##                       'measure_median': 'performance_median',
##                       'measure_CLs': 'performance_CLs'}
##
##ShuffledPerformancePlotSpecDict = {'measure': 'performance_median',
##                       'measure_median': 'performance_median',
##                       'measure_CLs': 'performance_CLs'}
##
##MutInfoPlotSpecDict = {'measure': 'mutual_info',
##                       'measure_median': 'mutual_info_median',
##                       'measure_CLs': 'mutual_info_CLs'}
##
##ShuffledMutInfoPlotSpecDict = {'measure': 'mutual_info_median',
##                       'measure_median': 'mutual_info_median',
##                       'measure_CLs': 'mutual_info_CLs'}
# 
## Plot performance dependence on increasing peri-event window span
## Define figure name
#drive, path_and_file = os.path.splitdrive(PathToBehavFile)
#path, file = os.path.split(path_and_file)
#FigureTitle = file[:-7]
# 
## Initialize figure
#fig1, axs1 = plt.subplots()
#fig1.suptitle(FigureTitle)
# 
## Plot performance and performance control plots
#PlotSpecDict = {'measure': 'performance',
#                'measure_median': 'performance_median',
#                'measure_CLs': 'performance_CLs',
#                'measure_SE': 'performance_SE',
#                'color':'blue'}
# 
#GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
#                     axs1, 'fw_sliding')
# 
#PlotSpecDict = {'measure': 'performance_median',
#                'measure_median': 'performance_median',
#                'measure_CLs': 'performance_CLs',
#                'measure_SE': 'performance_SE',
#                'color':'lightblue'}
# 
#GenerateConfIntsPlot(EventsShuffled, EventsShuffled, PlotSpecDict, 
#                     axs1, 'fw_sliding')
# 
#axs1.set_xbound(lower=ParamsDict['BoundaryWindow'][0], 
#                upper=ParamsDict['BoundaryWindow'][1])
#axs1.set_ybound(lower=0.4, upper=1.)
#
## Plot mutual information and mutual information control plots
#fig2, axs2 = plt.subplots()
#fig2.suptitle(FigureTitle)
# 
#PlotSpecDict = {'measure': 'mutual_info',
#                'measure_median': 'mutual_info_median',
#                'measure_CLs': 'mutual_info_CLs',
#                'measure_SE': 'mutual_info_SE',
#                'color':'blue'}
# 
#GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
#                     axs2, 'fw_sliding')
# 
#PlotSpecDict = {'measure': 'mutual_info_median',
#                'measure_median': 'mutual_info_median',
#                'measure_CLs': 'mutual_info_CLs',
#                'measure_SE': 'mutual_info_SE',
#                'color':'lightblue'}
# 
#GenerateConfIntsPlot(EventsShuffled, EventsShuffled, PlotSpecDict, 
#                     axs2, 'fw_sliding')
# 
#axs2.set_xbound(lower=ParamsDict['BoundaryWindow'][0], 
#                upper=ParamsDict['BoundaryWindow'][1])
#axs2.set_ybound(lower=0., upper=1.)


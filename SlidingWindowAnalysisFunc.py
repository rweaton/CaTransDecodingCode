#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:33:08 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd

import PeriEventTraceFuncLib as PETFL
import CalciumImagingFluorProcessing as CIFP

from collections import defaultdict

import os

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.messagebox import askyesno


# Define ParamsDict, the dictionary that contains the parameters for
# PLS decoding, Bootstrapping, Shuffle control processes.
# Peripheral target entry events

# Define a dictionary to contain keyed access to sliding window analysis parameters
ParamsDict = defaultdict(dict)

# List names of behavioral targets to decode.
# Right-arm-to-right-target, left-arm-to-left-target task
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']

# Right arm to both right and left targets
ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Left arm to both right and left targets
#ParamsDict['RefEventsList'] = ['M7T0_Entry_ts', 'M7T1_Entry_ts']

# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set peri-event extraction window (seconds)
#ParamsDict['BoundaryWindow'] = [-2., 2.]
ParamsDict['BoundaryWindow'] = [-2., 2.]

# Set width of sliding window (seconds)
ParamsDict['WindowWidth'] = 0.4

# Set increment size to slide window
ParamsDict['StepWidth'] = 0.1

# Specify the number of latent factors to include for the Partial 
# Least Squares classifier
ParamsDict['NumLatents'] = 5

# Specify the number of repetitions to perform for bootstrapping (Monte Carlo
# non-parametric and shuffled procedures)
ParamsDict['NumRepetitions'] = 30

# Set the confidence level for calculating confidence intervals
ParamsDict['ConfLevel'] = 0.95

# Specify a refractory window, relative to target entry, for removing target 
# re-entry occurences that followed shortly after the initial entry.  That is,
# specify the minimum time allowed between successive behavioral events (seconds).
# Comment out to include all behavior events.
#ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)

def SlidingWindowAnalysisFunc(BehavDict, CellFluorTraces_Frame, SavePath, ParamsDict):
    
    # Pack into a dict the two lists that specify the names and assigned 
    # numerical values of the behavioral event types listed in ParamsDict
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
        
    # Generate a set of windows used for subsequent "sliding window" decoding.
    ArrayOfSlidWindows = PETFL.SlidingWindowGen(ParamsDict['BoundaryWindow'], 
                                                ParamsDict['StepWidth'], 
                                                ParamsDict['WindowWidth'])
    
    ###### Behavioral dataframe generation/handling ######
    # If a path is given in in place of a dictionary for the second arg, then 
    # assemble the dictionary from the specfied filepath string.      
    if (type(BehavDict) == str):
        
        # Generate the (unfiltered) behavior dictionary comprised of  event 
        # occurence time stamps.
        
        # Set the string argument as the filepath 
        PathToBehavFile = BehavDict
        
        # Assemble the dictionary of event occurence timestamps
        BehavDict = PETFL.BehavDictGen(PathToBehavFile)
        
    # Otherwise, use the dictionary given in the argument.
    elif (type(BehavDict) == dict):
        
        pass
    
    # If the RelativeTolWindow field exists in ParamsDict, detect rapid repeats 
    # within each event list and then remove them
    if ParamsDict['RelativeTolWindow']:

        # Generate a filter for keeping only initial target entries, and 
        # removing rapid repeat entries that might have followed.
        EventFilters = PETFL.RemoveRepeatTargetEntries(BehavDict, 
                                            RefEventsDict['RefEventsList'], 
                                            ParamsDict['RelativeTolWindow'])
        # Remove repeat events
        for ef in EventFilters:
            
            BehavDict[ef] = BehavDict[ef][EventFilters[ef]] 
    
    ###### Fluorescence dataframe generation/handling ######
    # If a path is given in in place of a dataframe for the second arg, then 
    # import the dataframe from filepath string.  
    if type(CellFluorTraces_Frame) == str:
        
        # Generate the dataframe comprised of cell fluorescence traces.        
        
        # Set the string argument as the filepath 
        PathToFluorFile = CellFluorTraces_Frame
        
        # Assemble the dataframe of cell fluorescence traces.
        CellFluorTraces_Frame = PETFL.CellFluorTraces_FrameGen(PathToFluorFile)
            
        # Assemble dataframe of cell centroids.
        CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)
        
    # Otherwise, use the dataframe given in the argument.
    elif type(CellFluorTraces_Frame) == pd.DataFrame:
        
        pass
    
    # Remove events that occured either too early, or too late, with respect to 
    # the timeframe of the calcium fluorescence record. That is, events that 
    # would otherwise cause extraction of surrouding trace snippets that would
    # be incomplete (not of full length).
    BehavDict = PETFL.EarlyAndLateBehavEventRemoval(BehavDict, 
                                              CellFluorTraces_Frame, 
                                              RefEventsDict, 
                                              ParamsDict['BoundaryWindow'])
    
    # Initialize the array to contain the time domain boundaries of the 
    # series of sliding windows.
    (NumDomains, _) = ArrayOfSlidWindows.shape
    
    # Initialize an empty array to contain output dictionaries from the 
    # decoder cross-validation perfomance and Monte Carlo bootstrap routines
    Performance = np.empty((NumDomains,), dtype=dict)
    ConfInts = np.empty((NumDomains,), dtype=dict)
    EventsShuffled = np.empty((NumDomains,), dtype=dict)
    
    #######  Begin sliding window decoding ######
    # Iterate over windows list.  On each iteration perform PLS decoding on the
    # observed traces-set/target-type pair in the current window. Also perform 
    # Monte Carlo simulation in the domain to generate bootstrapped sampling
    # distributions.
    for i in np.arange(0, NumDomains):
     
        # Extract peri-event fluorescence within (event-relative) time domain
        # of current window.
        PeriEventExtractorDict = PETFL.PeriEventExtractor_Trace(BehavDict, 
                                        CellFluorTraces_Frame, RefEventsDict, 
                                        ArrayOfSlidWindows[i])
        
        # Write Peri-event activity element of dict to its own array
        #PEA_Array = PeriEventExtractorDict['PEA_Array']
        #(NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
        
        # Generate a set of indices to test the inclusion portion of the 
        # performance code.
        #InclusionSet = np.random.randint(0, high=NumTotalTrials, 
        #                                 size=(NumTotalTrials,))
    
        # Perform Partial Least Squares decoding using activity of current
        # (event-relative) time domain
        Performance[i] = PETFL.PLS_DecoderPerformance(PeriEventExtractorDict, 
                                                ParamsDict['NumLatents'])
        
        # Include the time domain of current window in the current output dict
        Performance[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
    
        # Run non-parametric bootstrap of activity to simulate sampling
        # distribution of PLS classifier performance.
        ConfInts[i] = PETFL.PLS_MonteCarlo(PeriEventExtractorDict,
                                     ParamsDict['NumLatents'], 
                                     ParamsDict['NumRepetitions'], 
                                     ParamsDict['ConfLevel'],
                                     Performance[i])
        
        # Include the time domain of current window in the current output dict        
        ConfInts[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})

        # Run PLS classifier after shuffling event type labels to simulate 
        # sampling distribution of chance associations
        EventsShuffled[i] = PETFL.PLS_Shuffle(PeriEventExtractorDict, 
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        
#        EventsShuffled[i] = PETFL.PLS_Shuffle2(PeriEventExtractorDict, 
#                                        ParamsDict['NumLatents'], 
#                                        ParamsDict['NumRepetitions'], 
#                                        ParamsDict['ConfLevel'])

        # Include the time domain of current window in the current output dict
        EventsShuffled[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
        
    #####################################
    ########   Start shelving   #########
    #####################################
    
    # Assemble save path
    # get root directory of save path from path to calcium data
    # SavePath is an arguement above.  The following script requires it 
    exec(open('./ShelveWorkspaceScript.py').read())
    
    ###########################################################################
    ######################        END OF FUNCTION       #######################
    ###########################################################################
    
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

CandidateSavePath = DefaultSavePath + os.path.sep + SaveFilename
result = askyesno(title="Save to default?", 
                  message="Save to default location:\n" + CandidateSavePath)

if result == True:
    
    # Check if default directory exists, if not, make it.
    if not os.path.exists(DefaultSavePath):
    
        os.makedirs(DefaultSavePath)
    
    # Set save path
    SavePath = CandidateSavePath

else:
    
    # Prompt user to enter the desired save path
    SavePath = asksaveasfilename(title='Set workspace save path',
                                 initialdir=DefaultSavePath,
                                 initialfile=SaveFilename)

# End file dialogs
root.withdraw()

# Run sliding window analysis using above-defined subroutine.
SlidingWindowAnalysisFunc(PathToBehavFile, PathToFluorFile, SavePath, ParamsDict)
#SlidingWindowAnalysisFunc(PathToBehavFile, FluorDataframe_Combined, SavePath, ParamsDict)
#SlidingWindowAnalysisFunc(BehavDict_Combined, CaImag_df_Combined, SavePath, ParamsDict)

###############################################################################
###################        END EXECUTION SCRIPT        ########################
###############################################################################

############ Miscellany from development evolution ############################

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


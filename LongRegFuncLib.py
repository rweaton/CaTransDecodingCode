#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 19:11:22 2020

@author: thugwithyoyo
"""
#from PeriEventTraceFuncLib import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import PeriEventTraceFuncLib as PETFL

import os
import shelve
from collections import defaultdict

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.messagebox import askyesno

#### Longitudinal Registration function library #### 

def BuildLongRegDataframeFromCSV(**kwargs):
    
    # Acquire path of workspace to load.
    # Set defualt parent directories
    DataRootDir = '/home/thugwithyoyo/CaTransDecoding/LongitudinalRegistrationLists'

    #if kwargs['path_to_long_cellIDs'] is not None:
    if np.isin('path_to_long_cellIDs', np.array(list(kwargs.keys()))):
        
        lr_path = kwargs['path_to_long_cellIDs']
        
    else:
        
        # Start file dialogs 
        root = tk.Tk()
        
        # Prompt user to navigate to, and select, the longitudinal cell identities file
        lr_path = askopenfilename(title='Select longitudinal cell identities file',
                                      filetypes=[("csv files", "*.csv")],
                                      initialdir=DataRootDir)

        # End file dialogs
        root.withdraw()
    
    # Asssemble longitudinal registration information into pandas datatrame
    lr_frame = pd.read_csv(lr_path)
    
    return lr_frame

def BuildSessionDatasets(**kwargs):
   
    # Acquire path of workspace to load.
    # Set defualt parent directories
    DataRootDir = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData'
        
    #if kwargs['path_to_behav_file'] is not None:
    if np.isin('path_to_behav_file', np.array(list(kwargs.keys()))):
        
        # If a path is included as an argument, set it to the path variable
        PathToBehavFile = kwargs['path_to_behav_file']

    else:
        
        # Start file dialogs 
        root = tk.Tk()
        
        # Prompt user to navigate to, and select, the session behavior file
        PathToBehavFile = askopenfilename(title='Select session behavior file',
                                          filetypes=[("json files","*.json"), 
                                                     ("csv files", "*.csv")],
                                          initialdir=DataRootDir)
        
        # End file dialogs
        root.withdraw()
        
    #if kwargs['path_to_fluor_file'] is not None:
    if np.isin('path_to_fluor_file', np.array(list(kwargs.keys()))):

        PathToFluorFile = kwargs['path_to_fluor_file']
        
    else:
        
        # Start file dialogs 
        root = tk.Tk()
        
        # Prompt user to navigate to, and select, the session imaging file
        PathToFluorFile = askopenfilename(title='Select session fluorescence record file',
                                          filetypes=[("json files","*.json"),
                                                     ("csv files", "*.csv")],
                                          initialdir=DataRootDir)
        # End file dialogs
        root.withdraw()
        
    # Assemble the dataframe of cell fluorescence traces.
    CellFluorTraces_Frame = PETFL.CellFluorTraces_FrameGen(PathToFluorFile)
    
    # Assemble the dictionary of event occurence timestamps
    BehavDict = PETFL.BehavDictGen(PathToBehavFile)

    # Determine parent directory and filename from complete path.
    Drive, Path_and_file = os.path.splitdrive(PathToBehavFile)
    Path, File = os.path.split(Path_and_file)
    
    # Extract session id from filename
    SessionID = File[0:19]
    
    return CellFluorTraces_Frame, BehavDict, SessionID


def ElementLocator(TargetsArray, TestArray):
    
    # This program finds the location of TargetsArray elements within the 
    # array TestArray, provided that they exist.  More specifically, for each
    # element in TargetsArray, the routine searches for a match in TestArray
    # and returns the location of that match (within TestArray) by its array 
    # index. The program assumes that elements of TestArray are unique; if 
    # multiple matches are found (1-to-n mapping) only the index of the first 
    # match is retained.
    
    # Count number of elements in TargetsArray
    (NumTargetElements,) = TargetsArray.shape
    
    # Count number of elements in  TestArray
    (NumTestElements,) = TestArray.shape
    
    # Generate sequential list of indices to count locations of elements in
    # TestArray
    TestArrayIndices = np.arange(0, NumTestElements)
    
    # Initialize the array to contain index locations to return
    TargetsMap = np.empty((NumTargetElements,), dtype=int)
    
    # Iterate over elements within TargetsArray
    for t in np.arange(0, NumTargetElements):
        
        # For each element, find matches within TestArray
        Filt = (TestArray == TargetsArray[t])
        
        # Write the index location of first match to the map output variable
        TargetsMap[t] = TestArrayIndices[Filt][0]
        
    return TargetsMap

def CommonCellsExtractor(lr_frame, Session1, Session2):
    
    # Generate filters to extract session subsets
    S1_filt = (lr_frame['local_cellset_file'] == Session1)
    S2_filt = (lr_frame['local_cellset_file'] == Session2)
    
    # Obtain cell lists by global index from the two sessions
    S1_global_cellset = lr_frame['global_cell_index'][S1_filt]
    S2_global_cellset = lr_frame['global_cell_index'][S2_filt]
    
    # Obtain cell lists by local cell name from  the two sessions
    S1_local_cellset = lr_frame['local_cell_name'][S1_filt]
    S2_local_cellset = lr_frame['local_cell_name'][S2_filt]
    
    # Determine particular cells that are present in both lists
    common_global_cellset = np.intersect1d(S1_global_cellset.values, 
                                           S2_global_cellset.values)
    
    # Sort the unique values of the result in ascending order
    common_global_cellset = np.sort(common_global_cellset)
    
    # Generate a map to cells, that are common to both sessions, within the
    # session 1 lists
    S1_cells_map = ElementLocator(common_global_cellset, 
                                    S1_global_cellset.values)
    
    # Generate a map to cells, that are common to both sessions, within the
    # session 2 lists
    S2_cells_map = ElementLocator(common_global_cellset,
                                    S2_global_cellset.values)
    
    return {
            'S1_cells_map' : S1_cells_map,
            'S2_cells_map' : S2_cells_map,
            'S1_global_cellset' : S1_global_cellset.values[S1_cells_map],
            'S2_global_cellset' : S2_global_cellset.values[S2_cells_map],
            'S1_local_cellset' : S1_local_cellset.values[S1_cells_map], 
            'S2_local_cellset' : S2_local_cellset.values[S2_cells_map]
            }

def TestOnTrainedPerformance(Performance_train,
                             PeriEventExtractorDict_train, 
                             PeriEventExtractorDict_test):
        
    #TrainingSetPerfMat[i, j] = Performance_train['performance']
    # Apply the train dataset to the current model trained above
    TrainingSetPerf = \
        PETFL.ApplyPLSModel(Performance_train['B'], 
                            Performance_train['B_0'], 
                            PeriEventExtractorDict_train)['performance']            
    
    
    # Apply the test dataset to the current model trained above
    TestingSetPerf = \
        PETFL.ApplyPLSModel(Performance_train['B'], 
                            Performance_train['B_0'], 
                            PeriEventExtractorDict_test)['performance']

    return {
            'TrainingSetPerf': TrainingSetPerf,
            'TestingSetPerf': TestingSetPerf,          
            }
    
def FilterBehavDict(ParamsDict, BehavDict, CellFluorTraces_Frame):
    
    # Assemble events reference dictionary 
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
    
    # Remove events that occured either too early, or too late, with 
    # respect to the timeframe of the calcium fluorescence record. 
    # That is, events that would otherwise cause extraction of 
    # surrounding trace snippets that would be incomplete (not of full 
    # length).
    BehavDict = PETFL.EarlyAndLateBehavEventRemoval(BehavDict, 
                                              CellFluorTraces_Frame, 
                                              RefEventsDict, 
                                              ParamsDict['BoundaryWindow'])
    
    if ParamsDict['RelativeTolWindow']:
        
        # Generate a filter for keeping only initial target entries, and 
        # removing rapid repeat entries that might have followed.
        EventFilters = PETFL.RemoveRepeatTargetEntries(BehavDict, 
                                            RefEventsDict['RefEventsList'], 
                                            ParamsDict['RelativeTolWindow'])

        # Remove repeat events in the training behavioral set
        for ef in EventFilters:
            
            BehavDict[ef] = BehavDict[ef][EventFilters[ef]]
        
    return BehavDict

    
def FilterCellFluorTracesFrames(CommonCellsDict, CellFluorTraces_Frame_train, 
                                CellFluorTraces_Frame_test):
        
    # Include only cells listed in common cells dictionary that are also
    # listed in the cell traces dataframe for training.
    Mask1 = np.isin(CommonCellsDict['S1_local_cellset'], 
                    CellFluorTraces_Frame_train.columns.values)       
    
    # Include only cells listed in common cells dictionary that are also
    # listed in the cell traces dataframe for testing.
    Mask2 = np.isin(CommonCellsDict['S2_local_cellset'], 
                    CellFluorTraces_Frame_test.columns.values)

    # Set intersect the two masks to obtain the final selection mask for
    # the cell fluorescence sets for both the training and testing sessions
    SelectFilt = (Mask1 & Mask2)
        
    # Select common cell traces from TRAINING set            
    CellFluorTraces_Frame_train = CellFluorTraces_Frame_train[
            np.hstack([np.array(['Timestamps']), 
                       CommonCellsDict['S1_local_cellset'][SelectFilt]])]
    
    # Select common cell traces from TESTING set
    CellFluorTraces_Frame_test = CellFluorTraces_Frame_test[
            np.hstack([np.array(['Timestamps']), 
                       CommonCellsDict['S2_local_cellset'][SelectFilt]])]
         
    return SelectFilt, CellFluorTraces_Frame_train, CellFluorTraces_Frame_test

def RunTestOnTrained(ParamsDict, CommonCellsDict, 
                     CellFluorTraces_Frame_train, BehavDict_train, 
                     CellFluorTraces_Frame_test, BehavDict_test):
    
    # Pack into a dict the two lists that specify the names and assigned 
    # numerical values of the behavioral event types listed in ParamsDict.
    # This dictionary will be used in later calls to the PLS decoder.
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
    
    # Unpack size variables for array initialization
    NumShuffles = ParamsDict['NumRepetitions']
    NumLatents = ParamsDict['NumLatents']
    
    # Initialize storage variables.
    ShufflesListIndices = np.arange(0, NumShuffles)
    OutcomeShufTrainingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
    OutcomeShufTestingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
    CellShufTrainingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
    CellShufTestingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
    
    
#    # Include only cells listed in common cells dictionary that are also
#    # listed in the cell traces dataframe for training.
#    Mask1 = np.isin(CommonCellsDict['S1_local_cellset'], 
#                    CellFluorTraces_Frame_train.columns.values)       
#    
#    # Include only cells listed in common cells dictionary that are also
#    # listed in the cell traces dataframe for testing.
#    Mask2 = np.isin(CommonCellsDict['S2_local_cellset'], 
#                    CellFluorTraces_Frame_test.columns.values)
#
#    # Set intersect the two masks to obtain the final selection mask for
#    # the cell fluorescence sets for both the training and testing sessions
#    SelectFilt = (Mask1 & Mask2)
#    
#    # Count the number of common cells to be used for this combination
#    # of training and testing sessions.
#    NumOfCommonCells = np.sum(SelectFilt)
#    
#    # Select common cell traces from TRAINING set            
#    CellFluorTraces_Frame_train = CellFluorTraces_Frame_train[
#            np.hstack([np.array(['Timestamps']), 
#                       CommonCellsDict['S1_local_cellset'][SelectFilt]])]
#    
#    # Select common cell traces from TESTING set
#    CellFluorTraces_Frame_test = CellFluorTraces_Frame_test[
#            np.hstack([np.array(['Timestamps']), 
#                       CommonCellsDict['S2_local_cellset'][SelectFilt]])]

    # Remove early and late behavioral events with respect to fluorescence 
    # records, as well as rapid repeat events (if indicated), from TRAINING
    # and TESTING behavioral dicts
    BehavDict_train = FilterBehavDict(ParamsDict, BehavDict_train, 
                                      CellFluorTraces_Frame_train)
    
    BehavDict_test = FilterBehavDict(ParamsDict, BehavDict_test, 
                                      CellFluorTraces_Frame_test)    

    # Prune cell fluorescence dataframes so that they each contain a universally
    # common set of cells.
    SelectFilt, CellFluorTraces_Frame_train, CellFluorTraces_Frame_test = \
        FilterCellFluorTracesFrames(CommonCellsDict, CellFluorTraces_Frame_train, 
                                CellFluorTraces_Frame_test)
        
    # Count the number of common cells to be used for this combination
    # of training and testing sessions.
    NumOfCommonCells = np.sum(SelectFilt)    
    
    # Run peri-event extraction and assemble the dictionary to be 
    # used for PLS-regression training.
    PeriEventExtractorDict_train = \
        PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                       CellFluorTraces_Frame_train, 
                                       RefEventsDict, 
                                       ParamsDict['BoundaryWindow'])
        
    # Run peri-event extraction and assemble the dictionary to be 
    # used to test the PLS-regression model above.
    PeriEventExtractorDict_test = \
        PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                       CellFluorTraces_Frame_test, 
                                       RefEventsDict, 
                                       ParamsDict['BoundaryWindow'])
            
    # Compute the PLS regression matrices from the TRAINING dataset.
    Performance_train = \
        PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                     NumLatents)
        
    #TrainingSetPerfMat[i, j] = Performance_train['performance']
    # Apply the train dataset to the current model trained above
#    TrainingSetPerf = \
#        PETFL.ApplyPLSModel(Performance_train['B'], 
#                            Performance_train['B_0'], 
#                            PeriEventExtractorDict_train)['performance']            
#    
#    
#    # Apply the test dataset to the current model trained above
#    TestingSetPerf = \
#        PETFL.ApplyPLSModel(Performance_train['B'], 
#                            Performance_train['B_0'], 
#                            PeriEventExtractorDict_test)['performance']
    
    TestOnTrainedPerfDict = \
        TestOnTrainedPerformance(Performance_train,
                                 PeriEventExtractorDict_train, 
                                 PeriEventExtractorDict_test)
    
    # Determine number of trials and total features from current peri-event
    # activity array.
    (nTotalTrials, nTotalFeatures) = \
        PeriEventExtractorDict_train['PEA_Array'].shape
            
    for k in ShufflesListIndices:
        
        #### TRAINING SET ####
        # Shuffle target outcomes of TRAINING set and retrain
        # Generate a map to shuffle target outcomes
        ShuffledIndices = np.arange(0, nTotalTrials) 
        np.random.shuffle(ShuffledIndices)
        
        # Shuffle target outcome correspondence to be passed to PLS trainer
        PeriEventExtractorDict_train['TargetsVec'] = \
            PeriEventExtractorDict_train['TargetsVec'][ShuffledIndices]
            
        # Train PLS model on shuffled outcomes data
        ShufOutcomesPerformance_train = \
            PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                         NumLatents)
            
#        # Run TESTING set in shuffled-outcome model and record performance
#        OutcomeShufTestingSetPerfVec[k] = \
#            PETFL.ApplyPLSModel(ShufOutcomesPerformance_train['B'], 
#                                ShufOutcomesPerformance_train['B_0'], 
#                                PeriEventExtractorDict_test)['performance']
#            
#        # Record performance in storage matrix
#        OutcomeShufTrainingSetPerfVec[k] = \
#            ShufOutcomesPerformance_train['performance']
        
        # Run TRANING and TESTING sets in shuffled-outcome model and record 
        # performance
        TestOnTrainedTrialShufPerfDict = \
            TestOnTrainedPerformance(ShufOutcomesPerformance_train,
                                     PeriEventExtractorDict_train, 
                                     PeriEventExtractorDict_test)
        
        # Record performance in TRAINING and TESTING storage matrices
        OutcomeShufTrainingSetPerfVec[k] = \
            TestOnTrainedTrialShufPerfDict['TrainingSetPerf']
            
        OutcomeShufTestingSetPerfVec[k] =  \
            TestOnTrainedTrialShufPerfDict['TestingSetPerf']
        
        # Shuffle TRAINING set cell order and run model on the result
        # Shuffle order of cell fluorescence traces in the dataframe
        CellFluorTraces_Frame_train_shuf = \
            CellFluorTraces_Frame_train[
                    np.hstack([np.array(['Timestamps']), 
                               np.random.choice(CommonCellsDict['S1_local_cellset'][SelectFilt], 
                                                size=NumOfCommonCells, 
                                                replace=False)])]
        
        # Generate test data structure with shuffled cell order
        # of the TRAINING set
        PeriEventExtractorDict_train_shuf = \
            PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                            CellFluorTraces_Frame_train_shuf, 
                            RefEventsDict, 
                            ParamsDict['BoundaryWindow'])

        # Run model on shuffled data structure and record performance
#        CellShufTrainingSetPerfVec[k] = \
#            PETFL.ApplyPLSModel(Performance_train['B'], 
#                                Performance_train['B_0'], 
#                                PeriEventExtractorDict_train_shuf)['performance']
        
        #### TESTING SET ####
        # Shuffle TESTING set cell order and run model on the result
        # Shuffle order of cell fluorescence traces in the dataframe
        CellFluorTraces_Frame_test_shuf = \
            CellFluorTraces_Frame_test[
                    np.hstack([np.array(['Timestamps']), 
                               np.random.choice(CommonCellsDict['S2_local_cellset'][SelectFilt], 
                                                size=NumOfCommonCells, 
                                                replace=False)])]
        
        # Generate test data structure with shuffled cell order
        # of the TESTING set
        PeriEventExtractorDict_test_shuf = \
            PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                            CellFluorTraces_Frame_test_shuf, 
                            RefEventsDict, 
                            ParamsDict['BoundaryWindow'])

        # Run model on shuffled data structure and record performance
#        CellShufTestingSetPerfVec[k] = \
#            PETFL.ApplyPLSModel(Performance_train['B'], 
#                                Performance_train['B_0'], 
#                                PeriEventExtractorDict_test_shuf)['performance']

        TestOnTrainedCellShufPerfDict = \
            TestOnTrainedPerformance(Performance_train,
                                     PeriEventExtractorDict_train_shuf, 
                                     PeriEventExtractorDict_test_shuf)
            
        CellShufTrainingSetPerfVec[k] = \
            TestOnTrainedCellShufPerfDict['TrainingSetPerf']
        
        CellShufTestingSetPerfVec[k] = \
            TestOnTrainedCellShufPerfDict['TestingSetPerf']

    return {
            'NumOfCommonCells': NumOfCommonCells,
            'TrainingSetCVPerf': Performance_train['performance'],            
            'TrainingSetPerf': TestOnTrainedPerfDict['TrainingSetPerf'],
            'TestingSetPerf': TestOnTrainedPerfDict['TestingSetPerf'],
            'OutcomeShufTrainingSetPerfVec': OutcomeShufTrainingSetPerfVec,
            'OutcomeShufTestingSetPerfVec': OutcomeShufTestingSetPerfVec,
            'CellShufTrainingSetPerfVec': CellShufTrainingSetPerfVec,
            'CellShufTestingSetPerfVec': CellShufTestingSetPerfVec            
            }
    
def SlidingWindowRunTestOnTrained(ParamsDict, CommonCellsDict, 
                     CellFluorTraces_Frame_train, BehavDict_train, 
                     CellFluorTraces_Frame_test, BehavDict_test):
    
    # Urgh... endless scope problems with shelve.
    #ParamsDict = GlobalParamsDict
    #ParamsDict['BoundaryWindow'] = [-.5, .1]
    #ParamsDict['BoundaryWindow'] = [-2., 2.]
    
    # Pack into a dict the two lists that specify the names and assigned 
    # numerical values of the behavioral event types listed in ParamsDict
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}

        
    # Generate a set of windows used for subsequent "sliding window" decoding.
    ArrayOfSlidWindows = PETFL.SlidingWindowGen(ParamsDict['BoundaryWindow'], 
                                                ParamsDict['StepWidth'], 
                                                ParamsDict['WindowWidth'])

#    print('Boundary window: ' + str(ParamsDict['BoundaryWindow']))
#    print('Windows to process: ' + str(ArrayOfSlidWindows))
    
    # Initialize the array to contain the time domain boundaries of the 
    # series of sliding windows.
    (NumDomains, _) = ArrayOfSlidWindows.shape
    
    # Unpack size variables for array initialization
    NumShuffles = ParamsDict['NumRepetitions']
    NumLatents = ParamsDict['NumLatents']
    ConfLevel = ParamsDict['ConfLevel']
    
    # Initialize storage variables.
    ShufflesListIndices = np.arange(0, NumShuffles)
    
    # Remove early and late behavioral events with respect to fluorescence 
    # records, as well as rapid repeat events (if indicated), from TRAINING
    # and TESTING behavioral dicts
    BehavDict_train = FilterBehavDict(ParamsDict, BehavDict_train, 
                                      CellFluorTraces_Frame_train)
    
    BehavDict_test = FilterBehavDict(ParamsDict, BehavDict_test, 
                                     CellFluorTraces_Frame_test)
    
    # Prune cell fluorescence dataframes so that they each contain a universally
    # common set of cells.
    SelectFilt, CellFluorTraces_Frame_train, CellFluorTraces_Frame_test = \
        FilterCellFluorTracesFrames(CommonCellsDict, CellFluorTraces_Frame_train, 
                                CellFluorTraces_Frame_test)
    # Count the number of common cells to be used for this combination
    # of training and testing sessions.
    NumOfCommonCells = np.sum(SelectFilt)  
    
    # Initialize an empty array to contain output dictionaries from the 
    # decoder cross-validation perfomance and Monte Carlo bootstrap routines
    DecodingStabilityDictArray = np.empty((NumDomains,), dtype=dict)
    
    #######  Begin sliding window decoding ######
    # Iterate over windows list.  On each iteration perform PLS regression
    # of the TRAINING set and apply the TESTING set to the resultant model.  On
    # each iteration also shuffle target outcomes and cell fluorescence order
    # to verify model validity.
    for t in np.arange(0, NumDomains): 
        
        #ParamsDict['BoundaryWindow'] = list(ArrayOfSlidWindows[t])
        
        # Have a heart; don't keep the poor user in suspense.
        print('Processing time domain: ' + str(ArrayOfSlidWindows[t]))
        
        # Initialize storage arrays for shuffling output    
        OutcomeShufTrainingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
        OutcomeShufTestingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
        CellShufTrainingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
        CellShufTestingSetPerfVec = np.empty_like(ShufflesListIndices, dtype=float)
        
        # Extract peri-event fluorescence within (event-relative) time domain
        # of current window.
        PeriEventExtractorDict_train = \
            PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                           CellFluorTraces_Frame_train, 
                                           RefEventsDict, 
                                           ArrayOfSlidWindows[t])

        
        # Run peri-event extraction and assemble the dictionary to be 
        # used to test the PLS-regression model.
        PeriEventExtractorDict_test = \
            PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                           CellFluorTraces_Frame_test, 
                                           RefEventsDict, 
                                           ArrayOfSlidWindows[t])
            
        # Compute the PLS regression matrices from the TRAINING dataset.
        Performance_train = \
            PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                         NumLatents)
        
        # Bootstrap the confidence band and SEs around performance of current window
        ConfInts_train = \
            PETFL.PLS_MonteCarlo(PeriEventExtractorDict_train, NumLatents,
                                 NumShuffles, ConfLevel, Performance_train)
                                 

        # Evaluate performance for the current TRAINING set and TEST set 
        # combination
        TestOnTrainedPerfDict = \
            TestOnTrainedPerformance(Performance_train,
                                     PeriEventExtractorDict_train, 
                                     PeriEventExtractorDict_test)
    
        # Determine number of trials and total features from current peri-event
        # activity array.
        (nTotalTrials, nTotalFeatures) = \
            PeriEventExtractorDict_train['PEA_Array'].shape
            
        ### Begin shuffling block    
        for k in ShufflesListIndices:
        
            #### TRAINING SET ####
            # Shuffle target outcomes of TRAINING set and retrain
            # Generate a map to shuffle target outcomes
            ShuffledIndices = np.arange(0, nTotalTrials) 
            np.random.shuffle(ShuffledIndices)
            
            # Shuffle target outcome correspondence to be passed to PLS trainer
            PeriEventExtractorDict_train['TargetsVec'] = \
                PeriEventExtractorDict_train['TargetsVec'][ShuffledIndices]
                
            # Train PLS model on shuffled outcomes data
            ShufOutcomesPerformance_train = \
                PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                             NumLatents)
            
            # Run TRANING and TESTING sets in shuffled-outcome model and record 
            # performance
            TestOnTrainedTrialShufPerfDict = \
                TestOnTrainedPerformance(ShufOutcomesPerformance_train,
                                         PeriEventExtractorDict_train, 
                                         PeriEventExtractorDict_test)
            
            # Record performance in TRAINING and TESTING storage matrices
            OutcomeShufTrainingSetPerfVec[k] = \
                TestOnTrainedTrialShufPerfDict['TrainingSetPerf']
                
            OutcomeShufTestingSetPerfVec[k] =  \
                TestOnTrainedTrialShufPerfDict['TestingSetPerf']
                
            # Shuffle TRAINING set cell order and run model on the result
            # Shuffle order of cell fluorescence traces in the dataframe
            CellFluorTraces_Frame_train_shuf = \
                CellFluorTraces_Frame_train[
                        np.hstack([np.array(['Timestamps']), 
                                   np.random.choice(CommonCellsDict['S1_local_cellset'][SelectFilt], 
                                                    size=NumOfCommonCells, 
                                                    replace=False)])]
            
            # Generate test data structure with shuffled cell order
            # of the TRAINING set
            PeriEventExtractorDict_train_shuf = \
                PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                CellFluorTraces_Frame_train_shuf, 
                                RefEventsDict, 
                                ArrayOfSlidWindows[t])
            
            #### TESTING SET ####
            # Shuffle TESTING set cell order and run model on the result
            # Shuffle order of cell fluorescence traces in the dataframe
            CellFluorTraces_Frame_test_shuf = \
                CellFluorTraces_Frame_test[
                        np.hstack([np.array(['Timestamps']), 
                                   np.random.choice(CommonCellsDict['S2_local_cellset'][SelectFilt], 
                                                    size=NumOfCommonCells, 
                                                    replace=False)])]
            
            # Generate test data structure with shuffled cell order
            # of the TESTING set
            PeriEventExtractorDict_test_shuf = \
                PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                CellFluorTraces_Frame_test_shuf, 
                                RefEventsDict, 
                                ArrayOfSlidWindows[t])
    
            TestOnTrainedCellShufPerfDict = \
                TestOnTrainedPerformance(Performance_train,
                                         PeriEventExtractorDict_train_shuf, 
                                         PeriEventExtractorDict_test_shuf)
                
            CellShufTrainingSetPerfVec[k] = \
                TestOnTrainedCellShufPerfDict['TrainingSetPerf']
            
            CellShufTestingSetPerfVec[k] = \
                TestOnTrainedCellShufPerfDict['TestingSetPerf']            
                
        DecodingStabilityDictArray[t] = {
            'PeriEventDomain': ArrayOfSlidWindows[t],
            'NumOfCommonCells': NumOfCommonCells,
            'TrainingSetCVPerf': Performance_train['performance'],
            'TrainingSetCVPerfSE': ConfInts_train['performance_SE'],
            'TrainingSetPerf': TestOnTrainedPerfDict['TrainingSetPerf'],
            'TestingSetPerf': TestOnTrainedPerfDict['TestingSetPerf'],
            'OutcomeShufTrainingSetPerfVec': OutcomeShufTrainingSetPerfVec,
            'OutcomeShufTestingSetPerfVec': OutcomeShufTestingSetPerfVec,
            'CellShufTrainingSetPerfVec': CellShufTrainingSetPerfVec,
            'CellShufTestingSetPerfVec': CellShufTestingSetPerfVec             
                }
            
    return DecodingStabilityDictArray

def OrganizeSlidingWindowOutputForPlotting(DecodingStabilityDictArray, ParamsDict, axs):
    
    # FYI: The argument here is output from SlidingWindowRunTestOnTrained 
    # routine above.
    (NumDicts,) = DecodingStabilityDictArray.shape
    
    # Initialize vectors to plot
    TimeVec = np.empty((NumDicts,), dtype=float)
    TrainingSetCVPerfVec = np.empty((NumDicts,), dtype=float)
    TrainingSetCVPerfSEMVec = np.empty((NumDicts,), dtype=float)    
    TrainingSetPerfVec = np.empty((NumDicts,), dtype=float)
    TestingSetPerfVec = np.empty((NumDicts,), dtype=float)
    OutcomeShufTrainingSetPerf_Mean = np.empty((NumDicts,), dtype=float)
    OutcomeShufTrainingSetPerf_SEM = np.empty((NumDicts,), dtype=float)
    OutcomeShufTestingSetPerf_Mean = np.empty((NumDicts,), dtype=float)
    OutcomeShufTestingSetPerf_SEM = np.empty((NumDicts,), dtype=float)
    CellShufTrainingSetPerf_Mean = np.empty((NumDicts,), dtype=float)
    CellShufTrainingSetPerf_SEM = np.empty((NumDicts,), dtype=float)
    CellShufTestingSetPerf_Mean = np.empty((NumDicts,), dtype=float)
    CellShufTestingSetPerf_SEM = np.empty((NumDicts,), dtype=float)
    
    # Count number of repetitions in shuffled vectors.
    (NumReps,) = DecodingStabilityDictArray[0]['OutcomeShufTrainingSetPerfVec'].shape

    # Iterate through dict array.  Populate plot vectors
    for t in np.arange(0, NumDicts):
        
        TimeVec[t] = (DecodingStabilityDictArray[t]['PeriEventDomain'][0] + 
                     (DecodingStabilityDictArray[t]['PeriEventDomain'][1] - 
                      DecodingStabilityDictArray[t]['PeriEventDomain'][0])/2.)
        
        TrainingSetCVPerfVec[t] = DecodingStabilityDictArray[t]['TrainingSetCVPerf']
        
        TrainingSetCVPerfSEMVec[t] = DecodingStabilityDictArray[t]['TrainingSetCVPerfSE']
        
        TrainingSetPerfVec[t] = DecodingStabilityDictArray[t]['TrainingSetPerf']
        
        TestingSetPerfVec[t] = DecodingStabilityDictArray[t]['TestingSetPerf']
        
        OutcomeShufTrainingSetPerf_Mean[t] = np.mean(
                DecodingStabilityDictArray[t]['OutcomeShufTrainingSetPerfVec'])
        
        OutcomeShufTrainingSetPerf_SEM[t] = np.std(
                DecodingStabilityDictArray[t]['OutcomeShufTrainingSetPerfVec'])/np.sqrt(NumReps)
        
        OutcomeShufTestingSetPerf_Mean[t] = np.mean(
                DecodingStabilityDictArray[t]['OutcomeShufTestingSetPerfVec'])
        
        OutcomeShufTestingSetPerf_SEM[t] = np.std(
                DecodingStabilityDictArray[t]['OutcomeShufTestingSetPerfVec'])/np.sqrt(NumReps)
        
        CellShufTrainingSetPerf_Mean[t] = np.mean(
                DecodingStabilityDictArray[t]['CellShufTrainingSetPerfVec'])
        
        CellShufTrainingSetPerf_SEM[t] = np.std(
                DecodingStabilityDictArray[t]['CellShufTrainingSetPerfVec'])/np.sqrt(NumReps)
        
        CellShufTestingSetPerf_Mean[t] = np.mean(
                DecodingStabilityDictArray[t]['CellShufTrainingSetPerfVec'])
        
        CellShufTestingSetPerf_SEM[t] = np.std(
                DecodingStabilityDictArray[t]['CellShufTrainingSetPerfVec'])/np.sqrt(NumReps)
        
    # Initialize figure and plotting axis
    #figure, axs = plt.subplots(nrows=1, ncols=1)
    
    # Set colors for each trace type
    ShuffledTrialsColor = 'lightblue'
    ShuffledCellsColor = 'gray'
    TrainingSetColor = 'red'
    TestingSetColor = 'purple'
    TrainingSetCVColor = 'orange'
   
    # To start, plot underlying standard error bands.
    axs.fill_between(TimeVec, TrainingSetCVPerfVec + 
                     TrainingSetCVPerfSEMVec,
                     TrainingSetCVPerfVec - 
                     TrainingSetCVPerfSEMVec,
                     color=TrainingSetCVColor,
                     alpha=0.5)
    
#    axs.fill_between(TimeVec, OutcomeShufTrainingSetPerf_Mean + 
#                     OutcomeShufTrainingSetPerf_SEM, 
#                     OutcomeShufTrainingSetPerf_Mean - 
#                     OutcomeShufTrainingSetPerf_SEM, 
#                     color=ShuffledTrialsColor,
#                     alpha=0.5)
    
    axs.fill_between(TimeVec, OutcomeShufTestingSetPerf_Mean + 
                     OutcomeShufTestingSetPerf_SEM, 
                     OutcomeShufTestingSetPerf_Mean - 
                     OutcomeShufTestingSetPerf_SEM, 
                     color=ShuffledTrialsColor,
                     alpha=0.5)
    
#    axs.fill_between(TimeVec, CellShufTrainingSetPerf_Mean + 
#                     CellShufTrainingSetPerf_SEM, 
#                     CellShufTrainingSetPerf_Mean - 
#                     CellShufTrainingSetPerf_SEM, 
#                     color=ShuffledCellsColor,
#                     alpha=0.5)    

    axs.fill_between(TimeVec, CellShufTestingSetPerf_Mean + 
                     CellShufTestingSetPerf_SEM, 
                     CellShufTestingSetPerf_Mean - 
                     CellShufTestingSetPerf_SEM, 
                     color=ShuffledCellsColor,
                     alpha=0.5)
    
    # Plot means
    axs.plot(TimeVec, TrainingSetCVPerfVec, '.-',
             color=TrainingSetCVColor, label='Training set cross-validation')    
    
#    axs.plot(TimeVec, TrainingSetPerfVec, '.-',
#             color=TrainingSetColor, label='Training set on trained model')
    
    axs.plot(TimeVec, TestingSetPerfVec, '.-',
             color=TestingSetColor, label='Testing set on trained model')
    
#    axs.plot(TimeVec, OutcomeShufTrainingSetPerf_Mean, '.-',
#             color=ShuffledTrialsColor, label='Training set shuffled outcomes')
    
    axs.plot(TimeVec, OutcomeShufTestingSetPerf_Mean, '.-',
             color=ShuffledTrialsColor, label='Testing set shuffled outcomes')
    
#    axs.plot(TimeVec, CellShufTrainingSetPerf_Mean, '.-',
#             color=ShuffledCellsColor, label='Traning set shuffled cells')
    
    axs.plot(TimeVec, CellShufTestingSetPerf_Mean, '.-', 
             color=ShuffledCellsColor, label='Testing set shuffled cells')
    

    # Label axes, and set axes bounds
    axs.set_xbound(lower=ParamsDict['BoundaryWindow'][0], 
                    upper=ParamsDict['BoundaryWindow'][1])
    axs.set_ybound(lower=0., upper=1.)
    
#    axs.set_xlabel('time relative to reach entry (s)')
#    axs.set_ylabel('performance (%)')
    axs.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axs.set_yticklabels(['0', '', '50', '', '100'])
    
    # Include legend
    #axs.legend()
    
def LoadDecodingStabilityWorkspace():

    DataRootDir = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/'
    root = tk.Tk()
    RestoreFilePath = askopenfilename(title='Select shelved workspace file',
                                          filetypes=[("shelve files","*.dat")],
                                          initialdir=DataRootDir)
    root.withdraw()

    # Determine parent directory and filename from complete path.
    Drive, Path_and_file = os.path.splitdrive(RestoreFilePath)
    Path, File = os.path.split(Path_and_file)

    try:
        
        exec(open('./RestoreShelvedWorkspaceScript.py').read())
    
    except:
        
        print('Error when restoring shelved workspace.  Will attempt to continue')

###### Quick plotter ############## 
def BatchPlotter(DecodingStabilityOutput, ParamsDict):
    
    ShapeTup = DecodingStabilityOutput.shape
    
    if (len(ShapeTup) == 1):

      DecodingStabilityOutput = np.expand_dims(DecodingStabilityOutput, axis=1)
        
    
    (NumRows, NumCols) = DecodingStabilityOutput.shape
    
    fig, axType = plt.subplots(NumRows, NumCols, figsize=(11,8))
    
    if type(axType) != np.ndarray:
        
        axs = np.empty((1,1), dtype=type(axType))
        axs[0,0] = axType
        
    else:
        
        axs = axType
    
    for i in np.arange(0, NumRows):
        
        for j in np.arange(0, NumCols):
                           
                OrganizeSlidingWindowOutputForPlotting(DecodingStabilityOutput[i,j], 
                                                       ParamsDict, axs[i, j])
                
                if i == (NumRows - 1):
                    
                    axs[i, j].set_xlabel('time relative to reach entry (s)')
                    
                if j == 0:
                    
                    axs[i, j].set_ylabel('performance (%)')
                    
                if ((i == (NumRows - 1)) and (j == (NumCols-1))):
                    
                    axs[i,j].legend()
                
    #plt.figure()

#### Debug script #####
LoadDecodingStabilityWorkspace()
i = 1
j = 0
SingleOutput = np.empty((1,), dtype=np.ndarray)
SingleOutput[0] = DecodingStabilityOutput[i, j]
BatchPlotter(SingleOutput, ParamsDict)

axs = plt.gca()
axs.set_title('Training set: ' + DecodingSessionsList[i] + ', Testing set: ' 
              + DecodingSessionsList[j])

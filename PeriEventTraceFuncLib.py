#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 22:07:25 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd
from collections import defaultdict

import CalciumImagingBehaviorAnalysis as CIBA
import CalciumImagingFluorProcessing as CIFP
import PartialLeastSquares as PLS
import hwfunclib_limited as HFL

import os
import shelve
from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt

#filename, file_extension = os.path.splitext('/path/to/somefile.ext')

# Generate dictionary of reference events
def BehavDictGen(PathToBehavFile):
    
    # Routine to take Plexon-generated behavioral data (in the form of event
    # timestamps) and parse them into a python dictionary usable by processing/
    # decoding routines in this library.  
    
    # Determine particular subroutine to parse behavioral events.
    _, file_extension = os.path.splitext(PathToBehavFile)
    
    # Examine file extension to determmine the corresponding parser to call.
    # Subroutines are contained in the script library: CalciumImagingBehaviorAnalysis.py
    if (file_extension == '.csv'):
        
        BehavDict = CIBA.CinePlexCSV_parser(PathToBehavFile)
        
    elif (file_extension == '.json'):
        
        BehavDict = CIBA.CinePlexFromMatJSON_parser(PathToBehavFile)
        
    else:
        
        # Report an error if the given file is not either as .csv or .json
        print('Behavioral events: unrecognized file format.')
        
    return BehavDict

# Generate dataframe of calcium events
def CellFluorTraces_FrameGen(PathToFluorFile):
    
    # Routine to take Calcium fluorescence data, written in the form of a json 
    # file, and parse it into a pandas dataframe. The dataframe format is required
    # for processing/decoding routines in this libarary to access calcium trace
    # data.
    
    # Determine particular subroutine to parse calcium traces.
    _, file_extension = os.path.splitext(PathToFluorFile)
    
    # Examine file extension to determine the corresponding parser to call.  
    # Subroutines are contained in the script library: CalciumImagingFluorProcessing.py
    if (file_extension == '.csv'):
        
        CellFluorTraces_Frame = CIFP.IDPS_CellFluorTracesParser(PathToFluorFile)

    elif (file_extension == '.json'):
        
        CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)
        
    else:
        
        # Report an error if the given file is not either as .csv or .json
        print('Calcium fluorescence traces: unrecognized file format.')

    return CellFluorTraces_Frame

def RemoveRepeatTargetEntries(BehavDict, RefEventsList, RelativeTolWindow):

    # This routine removes timestamps that follow earlier occurances within
    # a specified duration defined by the argument RelativeTolWindow.  Timestamps
    # are listed in BehavDict and the particular event timestamp lists to filter
    # are listed in RefEventsList.
    
    # Peripheral target entry events
    #RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

    # Count the number of lists to process.
    #BehavDict = BehavDictGen(PathToBehavFile)
    nEventListsToFilter = len(RefEventsList)
    
    # Initialize dict to contain filter lists
    #RapidRepeatDict = np.empty((nEventListsToFilter,), dtype=dict)
    EventFilters = defaultdict(dict)

    # Iterate over the collection of event types.  Call the routine EventComparator
    # to detect "rapid-repeat" timestamps and then RemoveEventsInTolWindow...
    # to remove them.  Write resulting filtered lists to EventFilters dict
    # initialized above.    
    for i in np.arange(0, nEventListsToFilter):
        
        #ComparatorDict = CIBA.EventComparator(tlist1, tlist2, tol_window)
        RapidRepeatDict =  CIBA.EventComparator(
                                        BehavDict[RefEventsList[i]].values, 
                                        BehavDict[RefEventsList[i]].values, 
                                        RelativeTolWindow)
        
        EventFilters[RefEventsList[i]] = CIBA.RemoveEventsInTolWindowFiltGen(
                                            RapidRepeatDict)

    return EventFilters

def EarlyAndLateBehavEventRemoval(BehavDict, CellFluorTraces_Frame, 
                                  RefEventsDict, BoundaryWindow):
    
    # This routines removes events that occured either too early, or too late, 
    # with respect to the timeframe of the calcium fluorescence record, for 
    # proper extraction of their surrounding snippets of trace activity. That 
    # is, events that would otherwise cause extraction of surrouding trace 
    # snippets that would be incomplete (not of full length).
    
    # Count the number of event types to be decoded.
    EventIndices = np.arange(0, len(RefEventsDict['RefEventsList']))
    
    # Calculate earlymost and latemost thesholds for behavioral events to have
    # full length snippets.   
    EarlyThresh = CellFluorTraces_Frame['Timestamps'].values[0] - BoundaryWindow[0]
    LateThresh = CellFluorTraces_Frame['Timestamps'].values[-1] - BoundaryWindow[1]
    
    # Iterate over collections of events of each type. For each collection,
    # generate a filter that keeps only events that occured after the early 
    # threshold (EarlyThresh) and a filter keeping events that occured before
    # the late threshold (LateThresh).  Write list changes to the behavior
    # dictionary.
    for i in EventIndices:
        
        # Generate a boolean filter for keeping within-threshold event occurances
        EventsFilt = (BehavDict[RefEventsDict['RefEventsList'][i]].values > EarlyThresh) & \
                     (BehavDict[RefEventsDict['RefEventsList'][i]].values < LateThresh)
        
        # Remove super-threshold event occurances from behavior dictionary
        BehavDict[RefEventsDict['RefEventsList'][i]] = \
            pd.Series(BehavDict[RefEventsDict['RefEventsList'][i]].values[EventsFilt])
        
    return BehavDict
    
def PeriEventExtractor_Trace(BehavDict, CellFluorTraces_Frame, RefEventsDict, BoundaryWindow):
    
    # This routine takes information from behavioral and imaging sources
    # and parses the into an input dict for the Partial Least Squares decoding
    # algorithm PLS_DecoderPerformance and other routines: PLS_MonteCarlo and
    # PLS_Shuffle (below).  Events listed in RefEventsDict contain the names
    # of behavioral events to be decoded and BoundaryWindow specifies the peri-
    # event domain around events that trace activity is to be extracted.
    
    # Count the number of event types to be decoded.
    EventIndices = np.arange(0, len(RefEventsDict['RefEventsList']))
    
#    # Determine particular subroutine to parse behavioral events.
#    _, file_extension = os.path.splitext(PathToBehavFile)
#    
#    if (file_extension == '.csv'):
#        
#        BehavDict = CIBA.CinePlexCSV_parser(PathToBehavFile)
#        
#    elif (file_extension == '.json'):
#        
#        BehavDict = CIBA.CinePlexFromMatJSON_parser(PathToBehavFile)
#        
#    else:
#        
#        print('Behavioral events: unrecognized file format.')
#        
#    # Determine particular subroutine to parse calcium traces.
#    _, file_extension = os.path.splitext(PathToFluorFile)
#    
#    if (file_extension == '.csv'):
#        
#        CellFluorTraces_Frame = CIFP.IDPS_CellFluorTracesParser(PathToFluorFile)
#
#    elif (file_extension == '.json'):
#        
#        CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)
#        
#    else:
#        
#        print('Calcium fluorescence traces: unrecognized file format.')
        
    # Initialize dict for storing peri-event activity; one key per reference event.
    PeriEventActivityArraysDict = defaultdict()
    
    # Initialize input argument arrays for Partial Least Squares
    PEA_Array = np.array([[]]).transpose()
    TargetsVec = np.array([[]], dtype=int).transpose()
    
    # Generate a dictionary of lists of indices for referencing trials for each
    # event.
    TrialIndicesByEventDict = defaultdict()
    nTrialsPerEvent = np.empty((EventIndices.size,), dtype=int)
    
    # Populate PLS input arrays by stacking vertically sets of peri-event activity
    # snippets for each event type.
    for i in EventIndices:
        
        # Extract peri-event snippets about each event of the current event type. 
        PeriEventActivityArraysDict[RefEventsDict['RefEventsList'][i]] = CIFP.PeriEventTraceExtractor(
                CellFluorTraces_Frame, BehavDict[RefEventsDict['RefEventsList'][i]].values, 
                BoundaryWindow)
        
        # Determine the number of events in the current set.
        (nTrialsPerEvent[i], nSamplesPerTrace) = \
            PeriEventActivityArraysDict[RefEventsDict['RefEventsList'][i]].shape
    
        # If this iteration is the first, then initialize empty arrays of appropriate
        # size for stacking.
        if i == 0:
            
            PEA_Array = np.empty((0, nSamplesPerTrace))
            TargetsVec = np.empty((0, 1))
        
        # Aquire number of total trials to generate index list for current set.
        (nTotalTrials, nSamplesPerTrace) = PEA_Array.shape    
        
        # Generate and store list of indices for current event set that indicate
        # trace row locations in the composite activity array PEA_Array. 
        TrialIndicesByEventDict[RefEventsDict['RefEventsList'][i]] = np.arange(nTotalTrials, 
                               nTotalTrials + nTrialsPerEvent[i])
        
        # Append current peri-activity set to the composite array.        
        PEA_Array = np.vstack((PEA_Array, 
                               PeriEventActivityArraysDict[RefEventsDict['RefEventsList'][i]]))
        
        # Each event type will be assigned its own unique integer.  Stack sets of 
        # event integer lables in the same manner as the PEA_Array above.
        RefEventVec = RefEventsDict['AssignedEventVals'][i]*np.ones(
                (nTrialsPerEvent[i], 1), dtype=int)
        
        TargetsVec = np.vstack((TargetsVec, RefEventVec))
    
    return {    
            'PEA_Array' : PEA_Array,
            'TargetsVec' : np.asarray(TargetsVec, dtype=int),
            'RefEventsDict' : RefEventsDict,
            'TrialIndicesByEventDict' : TrialIndicesByEventDict
            }

def PLS_DecoderPerformance(PeriEventExtractorDict, nLatents, **kwargs):
    
    # Routine that uses leave-one-out cross validation to compile a confusion
    # matrix of output from the the partial least squares decoding algorithm.
    # The routine takes as input the peri event activities and corresponding
    # target values within the dict PeriEventExtractorDict.  nLatents contains
    # the interger number of latent variables to form the Partial Least Squares
    # decoding model. The optional variable, trials_to_include, must contain
    # an array of integers each listing the corresponding index of a activity
    # target pair to be included in the PLS decoding model.  This option is
    # used when this routine is called during the performance/mutual infor-
    # mation bootstrap via the PLS_MonteCarlo routine below.
    
    # Unpack dictionary contents
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    TargetsVec = PeriEventExtractorDict['TargetsVec']
    RefEventsDict = PeriEventExtractorDict['RefEventsDict']
    TrialIndicesByEventDict = PeriEventExtractorDict['TrialIndicesByEventDict']
    
    # Index number of included events
    EventIndices = np.arange(0, len(RefEventsDict['RefEventsList']))
    
    # Initialize confusion matrix from which mutual information will be determined.
    ConfusionMatrix = np.zeros((len(RefEventsDict['RefEventsList']), 
                                len(RefEventsDict['RefEventsList'])), dtype=int)
        
    #### LEAVE ONE OUT DECODER PERFORMANCE EVALUATION
    # Iterate through each of the input traces in PEA_Array.  Select the current
    # trace to be the test trace and use the rest to compute the PLS model decoder
    # matrices (B, B_0).  Using the test trace make a prediction and record the 
    # outcome in the confusion matrix.
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    PEA_ArrayRowIndices = np.arange(0, nTotalTrials)
 
    # If the optional trials_to_include argument is present then set the trial
    # set list to its value, otherwise use all trials.
    if np.isin('trials_to_include', np.array(list(kwargs.keys()))):
        
        TrialSetIndices = kwargs['trials_to_include']
    
    else:
        
        TrialSetIndices = np.arange(0, nTotalTrials)
        
    # Iterate through each single trial in the peri-event activity array.  Find
    # the event type to which the current trial belongs and compute the regression
    # algorithm using the remaining (n - 1) set of trials.        
    for TestTrialIndex in TrialSetIndices:
    
        # Find the event pool that contains the current trial
        for Event in RefEventsDict['RefEventsList']:
            
            if np.isin(TestTrialIndex, TrialIndicesByEventDict[Event]):
            
                CurrentEvent = Event
                
        # Compute the PLS output matrices for the (remaining) complement set of 
        # activity traces
        # Acquire current test trace
        TestActivity = np.array([PEA_Array[TestTrialIndex, :]])
        
        # Generate filter that selects ALL OTHER traces from the complete set of traces
        ComplementFilt = ~(PEA_ArrayRowIndices == TestTrialIndex)
        
        # Run Partial Least Squares regression on the complement set of peri-event
        # activity traces.
        B, B_0, _ = PLS.PLS1(PEA_Array[ComplementFilt,:], TargetsVec[ComplementFilt,:], 
                          nLatents)
        
    
        # Using the resulting PLS model, predict the event during which the test
        # trace was recorded.
        ModelOutput = TestActivity @ B + B_0
        
        #### Start of the "Discriminant Analysis" portion of the procedure.
        # Initialize array to contain proximity values
        EuclideanDiffVals = np.empty((EventIndices.size,), dtype=float)
        
        # Iterate through events list calculating Euclidean distance on each
        # pass.  Write result to EuclideanDiffVals array.
        for i in EventIndices:
            
            EuclideanDiffVals[i] = HFL.EuclideanDistanceCalculator(np.array(
                    RefEventsDict['AssignedEventVals'][i]), ModelOutput)
        #EuclideanDiffVals = np.array(RefEventsDict['AssignedEventVals']) - ModelOutput
        
        # Obtain the index of the smallest distance. This is the target value
        # predicted by the PLS model.
        PredictedEventIndex = np.argmin(EuclideanDiffVals)
                
        # Obtain the value of the actual target.
        TrueEventIndex = EventIndices[RefEventsDict['RefEventsList'] == np.array(CurrentEvent)][0]
        
        # Increment the cooresponding element of the confusion matrix (e.g. correct
        # prediction (on diagonal), or incorrect prediction (off diagonal))
        ConfusionMatrix[TrueEventIndex, PredictedEventIndex] += 1
        
    ######## End of leave-one-out cross validation. #########
        
    # Compute the PLS regression matrices for the ENTIRE SET of (specified) trials.   
    B, B_0, _ = PLS.PLS1(PEA_Array[TrialSetIndices, :], 
                      TargetsVec[TrialSetIndices, :], nLatents)
    
    # Run the PLS decoder on the complete dataset.
    ModelOutput = PEA_Array[TrialSetIndices, :] @ B + B_0
    
    # Initialize array to contain model predictions.
    Predicted = np.empty((nTotalTrials,1), dtype=int)
    
    # Iterate over trials
    for j in np.arange(0, nTotalTrials):
        
        # Iterate over event types.
        for i in EventIndices:
            
            # Calculate Euclidean distance
            EuclideanDiffVals[i] = HFL.EuclideanDistanceCalculator(np.array(
                    RefEventsDict['AssignedEventVals'][i]), ModelOutput[j])    
        
        # Acquire index of closest event.  This is the model's prediction.
        PredictedEventIndex = np.argmin(EuclideanDiffVals)
        
        # Record the prediction.
        Predicted[j,:] = RefEventsDict['AssignedEventVals'][PredictedEventIndex]
    
    return {    
            'confusion_matrix' : ConfusionMatrix,
            'performance' : HFL.PerformanceCalculator(ConfusionMatrix),
            'mutual_info' : HFL.MutualInformationFromConfusionMatrix(ConfusionMatrix),
            'targets' : np.asarray(TargetsVec, dtype=int),
            'predicted' : Predicted,
            'B' : B,
            'B_0' : B_0
            }
    
    #nIncorrects = np.sum(np.abs(np.squeeze(np.round(Predicted)) - TargetsVec.transpose()[0]))

##### Routine to run novel dataset using existing model trained via PLS-regression
def ApplyPLSModel(B, B_0, PeriEventExtractorDict, **kwargs):
    
    # Unpack dictionary contents
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    TargetsVec = PeriEventExtractorDict['TargetsVec']
    RefEventsDict = PeriEventExtractorDict['RefEventsDict']
    TrialIndicesByEventDict = PeriEventExtractorDict['TrialIndicesByEventDict']
    
    # Index number of included events
    EventIndices = np.arange(0, len(RefEventsDict['RefEventsList']))
    
    # Initialize confusion matrix from which mutual information will be determined.
    ConfusionMatrix = np.zeros((len(RefEventsDict['RefEventsList']), 
                                len(RefEventsDict['RefEventsList'])), dtype=int)
        
    # Determine number of trials (rows) and number of "features" (columns) from
    # shape of peri-event activity array.
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    #PEA_ArrayRowIndices = np.arange(0, nTotalTrials)
    
    # If the optional trials_to_include argument is present then set the trial
    # set list to its value, otherwise use all trials.
    if np.isin('trials_to_include', np.array(list(kwargs.keys()))):
        
        TrialSetIndices = kwargs['trials_to_include']
    
    else:
        
        TrialSetIndices = np.arange(0, nTotalTrials)
    
    # Determine the number of trials to be included    
    (nTrialSetIndices,) = np.array(TrialSetIndices).shape
    
    # Run the PLS decoder on the specified dataset.
    ModelOutput = PEA_Array[TrialSetIndices, :] @ B + B_0
    
    # Initialize array to contain proximity values
    EuclideanDiffVals = np.empty((EventIndices.size,), dtype=float)
        
    # Initialize array to contain model predictions.
    Predicted = np.empty((nTrialSetIndices, 1), dtype=int)
    
    # Initialize variable that specifies the actual event type on each trial
    # iteration.
    CurrentEvent = ''
    
    # Iterate over specified trials
    for j in np.arange(0, nTrialSetIndices):
            
        # Find the event pool that contains the current trial
        for Event in RefEventsDict['RefEventsList']:
            
            if np.isin(TrialSetIndices[j], TrialIndicesByEventDict[Event]):
            
                CurrentEvent = Event
                
        # Iterate over event types.
        for i in EventIndices:
            
            # Calculate Euclidean distance
            EuclideanDiffVals[i] = HFL.EuclideanDistanceCalculator(np.array(
                    RefEventsDict['AssignedEventVals'][i]), ModelOutput[j])
            
        # Acquire index of closest event.  This is the model's prediction.
        PredictedEventIndex = np.argmin(EuclideanDiffVals)
        
        # Obtain the value of the actual target.
        TrueEventIndex = EventIndices[RefEventsDict['RefEventsList'] == np.array(CurrentEvent)][0]
        
        # Increment the cooresponding element of the confusion matrix (e.g. correct
        # prediction (on diagonal), or incorrect prediction (off diagonal))
        ConfusionMatrix[TrueEventIndex, PredictedEventIndex] += 1        
        
        # Record the prediction.
        Predicted[j,:] = RefEventsDict['AssignedEventVals'][PredictedEventIndex]
        
        # Reset the CurrentEvent string for next trial iteration  
        CurrentEvent = '' 
    
    return {    
            'confusion_matrix' : ConfusionMatrix,
            'performance' : HFL.PerformanceCalculator(ConfusionMatrix),
            'mutual_info' : HFL.MutualInformationFromConfusionMatrix(ConfusionMatrix),
            'targets' : np.asarray(TargetsVec, dtype=int),
            'predicted' : Predicted
            }
    
##### PLS_MonteCarlo
def PLS_MonteCarlo(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel, ObservedPerfDict):
    
    # This function performs a non-parametric bootstrap of the partial least
    # squares decoding algorithm.  Input data (calcium fluorescence traces) and
    # decoding target values (target types) are contained in PeriEventExtractorDict--
    # output dict of PeriEventExtractor_Trace routine above. 
    
    # Calculate alphas and boundary indices from arguments
    alphas = np.empty((2,), dtype=int)
    alphas[0] = np.ceil(nRepetitions*(1. - ConfLevel)/2.)
    alphas[1] = np.floor(nRepetitions*(1. + ConfLevel)/2.)
    
    # Unpack needed dictionary contents
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    RefEventsDict = PeriEventExtractorDict['RefEventsDict']

    # Initialize arrays for storing simulated output distributions    
    nRefEvents = len(RefEventsDict['RefEventsList'])
    ConfusionMatrixSimDist = np.empty((nRefEvents, nRefEvents, nRepetitions),
                                      dtype=int)
    PerformanceSimDist = np.empty((nRepetitions,), dtype=float)
    MutualInfoSimDist = np.empty((nRepetitions,), dtype=float)
    
    # Generate set of indices from which to extract sets for each bootstrap 
    # iteration
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    #MasterIndicesSet = np.arange(0, nTotalTrials)
    
    # Repeat PLS decoding over the specified number of repetitions.
    for i in np.arange(0, nRepetitions):
        
        # Draw, with replacement, N random trials from the dataset (with N 
        # being equal to the total number of trials in the dataset)
        InclusionSet = np.random.randint(0, high=nTotalTrials, 
                                         size=(nTotalTrials,))
        
        # Invoking the optional keyword argument, trials_to_include, run the
        # PLS decoder on the bootstrapped-sampled dataset.
        PerformanceStruct = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    nLatents, trials_to_include=InclusionSet)
        
        # Write current confusion matrix to the collection of bootstrapped 
        # confusion matrices.
        ConfusionMatrixSimDist[:,:,i] = PerformanceStruct['confusion_matrix']
        
        # Write the current Performance dict to the collection of bootstrapped
        # performance dicts.
        PerformanceSimDist[i] = PerformanceStruct['performance']
        
        # Write the mutual information dict to the collection of bootstrapped
        # mutual informations dicts.
        MutualInfoSimDist[i] = PerformanceStruct['mutual_info']
        
    
    ### Compute statistics from bootstrapping. ###
    # Compute median values of the collection of confusion matrices.
    confusion_matrix_median = np.median(ConfusionMatrixSimDist, axis=-1)
    
    # Compute ConfLevel confidence interval limits from collection of bootstrapped
    # confidence matrices.
    confusion_matrix_CLs = np.sort(2*ObservedPerfDict['confusion_matrix'] - 
                                   np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas],
                                   axis=-1)
    
    # Compute ConfLevel confidence interval limits from collection of bootstrapped
    # of performance values.
    performance_CLs = np.sort(2.*ObservedPerfDict['performance'] - 
                              np.sort(PerformanceSimDist, axis=-1)[alphas], 
                              axis=-1)

    # Compute ConfLevel confidence interval limits from collection of bootstrapped
    # of mutual information values.    
    mutual_info_CLs = np.sort(2.*ObservedPerfDict['mutual_info'] - 
                              np.sort(MutualInfoSimDist, axis=-1)[alphas], 
                              axis=-1)
    
    # Compute the median of the collection of bootstrapped performance values.
    performance_median = np.median(PerformanceSimDist, axis=-1)
    
    # Compute the mean of the collection of bootstrapped performance values.
    performance_mean = np.mean(PerformanceSimDist, axis=-1)
    
    # Sort the list of bootstrapped performance values for bias visualization.
    performance_SortedSimDist = np.sort(PerformanceSimDist, axis=-1)
    
    # Obtain values of ConfLevel quantiles for bias visualization.
    performance_quantiles = performance_SortedSimDist[alphas]
    
    # Compute ConfLevel confidence interval limits about bootstrapped median
#    performance_CLs = np.sort(2.*performance_median - 
#                              performance_SortedSimDist[alphas], axis=-1)
    
    # Compute standard error of performance as the standard deviation of the 
    # performance bootstrap distribution.
    performance_SE = np.std(PerformanceSimDist, axis=-1, ddof=1)
    
    # Compute the boundaries of the standard error bars about the performance
    # median
#    performance_SE_bounds = np.array([performance_median - performance_SE,
#                                      performance_median + performance_SE])
    # Compute the boundaries of the standared error bars about the observed performance.
    performance_SE_bounds = np.array([ObservedPerfDict['performance'] - performance_SE,
                                      ObservedPerfDict['performance'] + performance_SE])
    
    # Compute the median value of the bootstrap distribution of mutual info.
    mutual_info_median = np.median(MutualInfoSimDist, axis=-1)
    
    # Compute the standard error of the mutual information as the standard deviation
    # of the bootstrap distribution.
    mutual_info_SE = np.std(MutualInfoSimDist, axis=-1, ddof=1)
    
    # Compute the boundaries of standard error bars about the observed mutual information
    # values.
    mutual_info_SE_bounds = np.array([ObservedPerfDict['mutual_info'] - mutual_info_SE,
                                      ObservedPerfDict['mutual_info'] + mutual_info_SE])
    
    # Return bootstrap statistics in dictionary form
    return {
            'confusion_matrix_median' : confusion_matrix_median,
            'confusion_matrix_CLs' : confusion_matrix_CLs,
            'performance_median' : performance_median,
            'performance_mean' : performance_mean,
            'performance_quanitles' : performance_quantiles,
            'performance_CLs' : performance_CLs,
            'performance_SE' : performance_SE,
            'performance_SortedSimDist' : performance_SortedSimDist,
            'performance_SE_bounds' : performance_SE_bounds,
            'mutual_info_median' : mutual_info_median,
            'mutual_info_CLs' : mutual_info_CLs,
            'mutual_info_SE' : mutual_info_SE,
            'mutual_info_SE_bounds' : mutual_info_SE_bounds   
            }
    
##### SHUFFLE CONTROL
# To test validity of decoding, shuffle event lables prior to running decoder
def PLS_Shuffle(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel):
    
    # This function is out-dated as it basically repeats the Monte Carlo bootstrap
    # process that is defined in PLS_MonteCarlo above.  This early version
    # of the routine is kept here for possible reference.  Strongly suggest 
    # using PLS_Shuffle2 instead.
    
    # Calculate alphas and boundary indices from arguments
    alphas = np.empty((2,), dtype=int)
    alphas[0] = np.ceil(nRepetitions*(1. - ConfLevel)/2.)
    alphas[1] = np.floor(nRepetitions*(1. + ConfLevel)/2.)
    
    # Unpack needed dictionary contents
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    RefEventsDict = PeriEventExtractorDict['RefEventsDict']

    # Initialize arrays for storing simulated output distributions    
    nRefEvents = len(RefEventsDict['RefEventsList'])
    ConfusionMatrixSimDist = np.empty((nRefEvents, nRefEvents, nRepetitions),
                                      dtype=int)
    PerformanceSimDist = np.empty((nRepetitions,), dtype=float)
    MutualInfoSimDist = np.empty((nRepetitions,), dtype=float)
    
    # Generate set of indices from which to extract sets for each bootstrap 
    # iteration
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    ShuffledIndices = np.arange(0, nTotalTrials)
    
    for i in np.arange(0, nRepetitions):
        
        np.random.shuffle(ShuffledIndices)
        PeriEventExtractorDict['TargetsVec'] = \
            PeriEventExtractorDict['TargetsVec'][ShuffledIndices]
        
        PerformanceStruct = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    nLatents)
        
        ConfusionMatrixSimDist[:,:,i] = PerformanceStruct['confusion_matrix']
        PerformanceSimDist[i] = PerformanceStruct['performance']
        MutualInfoSimDist[i] = PerformanceStruct['mutual_info']
            
    
#    return {
#            'confusion_matrix_median' : np.median(ConfusionMatrixSimDist, axis=-1),
#            'confusion_matrix_CLs' : np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas],
#            'performance_median' : np.median(PerformanceSimDist, axis=-1),
#            'performance_CLs' : np.sort(PerformanceSimDist, axis=-1)[alphas],
#            'performance_SE' : np.std(PerformanceSimDist, axis=-1, ddof=1),            
#            'mutual_info_median' : np.median(MutualInfoSimDist, axis=-1),
#            'mutual_info_CLs' : np.sort(MutualInfoSimDist, axis=-1)[alphas],
#            'mutual_info_SE' : np.std(MutualInfoSimDist, axis=-1, ddof=1)
#            }

    # Compute statistics from bootstrapping.
    confusion_matrix_median = np.median(ConfusionMatrixSimDist, axis=-1)
    confusion_matrix_CLs = np.sort(2.*confusion_matrix_median - 
                                   np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas], axis=-1)

    performance_median = np.median(PerformanceSimDist, axis=-1)
    performance_mean = np.mean(PerformanceSimDist, axis=-1)
    performance_SortedSimDist = np.sort(PerformanceSimDist, axis=-1)
    performance_quantiles = performance_SortedSimDist[alphas]
    performance_CLs = np.sort(2.*performance_median - 
                              performance_SortedSimDist[alphas], axis=-1)
    performance_SE = np.std(PerformanceSimDist, axis=-1, ddof=1)
    performance_SE_bounds = np.array([performance_median - performance_SE,
                                      performance_median + performance_SE])
    
    mutual_info_median = np.median(MutualInfoSimDist, axis=-1)
    mutual_info_CLs = np.sort(2.*mutual_info_median - 
                              np.sort(MutualInfoSimDist, axis=-1)[alphas], axis=-1)
    mutual_info_SE = np.std(MutualInfoSimDist, axis=-1, ddof=1)
    mutual_info_SE_bounds = np.array([mutual_info_median - mutual_info_SE,
                                      mutual_info_median + mutual_info_SE])
    
    return {
            'confusion_matrix_median' : confusion_matrix_median,
            'confusion_matrix_CLs' : confusion_matrix_CLs,
            'performance_median' : performance_median,
            'performance_mean' : performance_mean,
            'performance_quanitles' : performance_quantiles,
            'performance_CLs' : performance_CLs,
            'performance_SE' : performance_SE,
            'performance_SortedSimDist' : performance_SortedSimDist,
            'performance_SE_bounds' : performance_SE_bounds,
            'mutual_info_median' : mutual_info_median,
            'mutual_info_CLs' : mutual_info_CLs,
            'mutual_info_SE' : mutual_info_SE,
            'mutual_info_SE_bounds' : mutual_info_SE_bounds            
            }    


def PLS_Shuffle2(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel):
    
    # This version of partial least squares 
    # Calculate alphas and boundary indices from arguments
    alphas = np.empty((2,), dtype=int)
    alphas[0] = np.ceil(nRepetitions*(1. - ConfLevel)/2.)
    alphas[1] = np.floor(nRepetitions*(1. + ConfLevel)/2.)
    
    # Unpack needed dictionary contents
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    RefEventsDict = PeriEventExtractorDict['RefEventsDict']

    # Initialize arrays for storing simulated output distributions    
    nRefEvents = len(RefEventsDict['RefEventsList'])
    ConfusionMatrixSimDist = np.empty((nRefEvents, nRefEvents, nRepetitions),
                                      dtype=int)
    PerformanceSimDist = np.empty((nRepetitions,), dtype=float)
    MutualInfoSimDist = np.empty((nRepetitions,), dtype=float)
    
    # Generate set of indices from which to extract sets for each bootstrap 
    # iteration
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    ShuffledIndices = np.arange(0, nTotalTrials)
    
    # Shuffle order of referencing indices to destroy correspondence between
    # peri-event fluorescence activity snippits and the particular events 
    # around which the snippets were extracted.   
    np.random.shuffle(ShuffledIndices)
    PeriEventExtractorDict['TargetsVec'] = \
        PeriEventExtractorDict['TargetsVec'][ShuffledIndices]
    
    # Run PLS decoding on the shuffled dataset.  In this context, the output
    # performance we are calling "observed" for the MonteCarlo bootstrapping
    # to follow.
    ObservedPerfDict = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                nLatents)
    
    # Perform a non-paramentric bootstrap on the shuffled dataset using the 
    # PLS_MonteCarlo routine above.
    return PLS_MonteCarlo(PeriEventExtractorDict, nLatents, nRepetitions, 
                          ConfLevel, ObservedPerfDict)
        
    
#TargetsVec = TargetsVec.transpose()[0]
#np.random.shuffle(TargetsVec)
#TargetsVec = np.array([TargetsVec]).transpose()

#### The following function does not work as intended!!! ####
# See: ShelveWorkspaceScript.py 
#def ShelveWorkspace(SavePath):
#    
#    ScriptDir = os.getcwd()
#    
#    # Attempt to use GUI navigator dialog to guide save.  Drop for batch processing.
##    root = Tk()
##    root.filename =  filedialog.asksaveasfilename(
##                        initialdir = "/home/thugwithyoyo/Desktop",
##                        title = "Enter name of file to save", 
##                        filetypes = (("out files","*.out"), ("all files","*.*")))    
##    
##    if root.filename is None: # asksaveasfile return `None` if dialog closed with "cancel".
##        return
##    
##    LoadFilePath = root.filename
##    root.destroy()
#    
#    # Use path from function argument to change to save directory and write file
#    # with name contained in Filename.
#    (PathToFile, Filename) = os.path.split(SavePath)
#    
#    os.chdir(PathToFile)
#     
#    # Open a shelf object to contain workspace variables.
#    my_shelf = shelve.open(Filename)    
#    #PathToSaveFile = root.filename
#    #my_shelf = shelve.open(PathToSaveFile, 'n') # 'n' for new
#
#    # Iterate through list of global variables and write to file.
#    for key in dir():
#        
#        try:
#            
#            my_shelf[key] = globals()[key]
#            
#        except TypeError:
#            #
#            # __builtins__, my_shelf, and imported modules can not be shelved.
#            #
#            print('ERROR shelving: {0}'.format(key))
#            
#    #root.destroy()
#    
#    # Close shelf object after variables have been written to file
#    my_shelf.close()
#    
#    # Return to directory that contains this script.
#    os.chdir(ScriptDir)

#### The following function does not work as intended!!! ####
# See: RestoreShelvedWorkspaceScript.py
#def RestoreShelvedWorkspace(RestoreFilePath):
#    
#    ScriptDir = os.getcwd()
#    
##    root = Tk()
##    root.filename =  filedialog.askopenfilename(
##            initialdir = "/home/thugwithyoyo/Desktop",
##            title = "Select file to open",
##            filetypes = (("dat files","*.dat"),("all files","*.*")))
##
##    RestoreFilePath = root.filename
##    root.destroy()
#    (PathToFile, Filename) = os.path.split(RestoreFilePath)
#    
#    os.chdir(PathToFile)
#     
#    my_shelf = shelve.open(os.path.splitext(Filename)[0])
#    
#    for key in my_shelf:
#        
#        globals()[key]=my_shelf[key]
#        
#    my_shelf.close()
#    
#    os.chdir(ScriptDir)
    
def GenerateConfIntsPlot(ArrayOfConfIntsDicts, ArrayOfPerformanceDicts, PlotSpecDict, AxesHandle, Direction):
    
    # Routine to generate plots of decoding performance or mutual information
    # values over (event-relative) time.  The routine takes as arguments
    # output arrays from either SlidingWindowAnalysisFunc or 
    # VariableWindowAnalysisFunc routines.  Generated plots contain a solid line
    # to depict the observed results (performance or mutual information) and
    # corresponding error-values (either as confidence intervals or standard 
    # errors) as an underlying filled band.  The input argument PlotSpecDict 
    # contains directives about which items in the data arrays (ArrayofConfIntsDicts
    # and ArrayOfPerformanceDicts) are to used to build the plot. AxesHandle contains
    # a handle to the pre-exisiting axes object on which the the traces will be drawn.
    
    # Count number of datapoints and confidence intervals to plot.
    (nDicts,) = ArrayOfConfIntsDicts.shape
    
    # Determine the spacing between successive datapoints.
    StepSize = np.max(np.abs(ArrayOfConfIntsDicts[1]['PeriEventDomain'] - 
                             ArrayOfConfIntsDicts[0]['PeriEventDomain']))
    
    # Initialize data summary arrays to send to plotting routines.
    X = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)
    Y_median = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)
#    CI_Bars = np.empty((nDicts,2), dtype=float)
    Y_BarBounds = np.empty((nDicts,2), dtype=float)
    Y_test = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)

    # Obtain the y-axis name. If "_median" substring is at the end of the 
    # measure name, remove it.  Obtain the name of the trace for legend.
    OrdinateName = PlotSpecDict['measure']
    
    if OrdinateName.endswith('_median'):
        OrdinateName = OrdinateName[:-7]
        TraceType = 'Shuffled outcomes'
    else:
        TraceType = 'Observed outcomes'
    
    # Generate plot of performance measure over event-relative time plot with error band.
    # In this locations of datapoints along the time axis are defined and labled
    # for time axis and plot title.
    for i in np.arange(0, nDicts):
        
        # Use this clause if domain is backwards-cumulative.
        if (Direction == 'bw'):
            
            X[i] = np.min(ArrayOfConfIntsDicts[i]['PeriEventDomain']) + 0.5*StepSize
            TraceLabel = 'neg. cumulative span from '+ \
                str(ArrayOfConfIntsDicts[i]['PeriEventDomain'][1]) +' sec.' \
                + '('+ TraceType + ')'
            TitleLabel = ' dependence on span of activity window'
        
        # Use this clause if domain is forwards-cummulative.
        elif (Direction == 'fw'):
            
            X[i] = np.max(ArrayOfConfIntsDicts[i]['PeriEventDomain']) - 0.5*StepSize
            TraceLabel = 'pos. cumulative span from '+ \
                str(ArrayOfConfIntsDicts[i]['PeriEventDomain'][0]) +' sec.' \
                + '('+ TraceType + ')'                
            TitleLabel = ' dependence on span of activity window'

        # Use this clause if domain is a foward-sliding window.
        elif (Direction == 'fw_sliding'):
            
            X[i] = np.mean(np.array(ArrayOfConfIntsDicts[i]['PeriEventDomain']))
            WindowWidth = np.diff(np.array(ArrayOfConfIntsDicts[0]['PeriEventDomain']))[0]
            TraceLabel = TraceType
            TitleLabel = ': window of width '+str(round(WindowWidth, 3))+' sec. slid over ' \
                            +str(round(StepSize, 3))+' sec. increments.'
            
#        Y_median[i] = ArrayOfConfIntsDicts[i]['performance_median']
        Y_median[i] = ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']]
        
#        CI_Bars[i,:] = np.abs(ArrayOfConfIntsDicts[i]['performance_CLs'] - 
#                              ArrayOfConfIntsDicts[i]['performance_median'])
#        CI_Bars[i,:] = np.abs(ArrayOfConfIntsDicts[i][PlotSpecDict['measure_CLs']] - 
#                              ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']])
        
        # Plot Standard Error (standard deviation of sampling distribution of 
        # the statistic) rather than Confidence Intervals
#        CI_Bars[i,:] = np.abs(ArrayOfConfIntsDicts[i][PlotSpecDict['measure_bars']] - 
#                              ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']])

#        CI_Bounds[i,:] = ArrayOfConfIntsDicts[i]['performance_CLs']
        Y_BarBounds[i,:] = ArrayOfConfIntsDicts[i][PlotSpecDict['measure_bars']]        
        
        # Plot band centered around median of sampling distribution
#        CI_Bounds[i,0] = (ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']] - 
#                         ArrayOfConfIntsDicts[i][PlotSpecDict['measure_SE']])
        
#        CI_Bounds[i,1] = (ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']] + 
#                         ArrayOfConfIntsDicts[i][PlotSpecDict['measure_SE']])
        
        # Plot band centered around observed performance values
#        CI_Bounds[i,0] = (ArrayOfPerformanceDicts[i][PlotSpecDict['measure']] - 
#                         ArrayOfConfIntsDicts[i][PlotSpecDict['measure_bars']])
#        
#        CI_Bounds[i,1] = (ArrayOfPerformanceDicts[i][PlotSpecDict['measure']] + 
#                         ArrayOfConfIntsDicts[i][PlotSpecDict['measure_bars']])        
#        Y_test[i] = ArrayOfPerformanceDicts[i]['performance']
        Y_test[i] = ArrayOfPerformanceDicts[i][PlotSpecDict['measure']]
        
    #fig = plt.figure()
    
    # Plot error band on specified axes object.
    #AxesHandle.errorbar(X, Y_median, yerr=CI_Bars.transpose(), lolims=True, uplims=True, label=TraceLabel)
    AxesHandle.fill_between(X, Y_BarBounds.transpose()[0,:], Y_BarBounds.transpose()[1,:], 
                            label=TraceLabel, alpha=0.7, color=PlotSpecDict['color'])
    
    # Plot observed performance measure on specified axes object.
    AxesHandle.plot(X, Y_test, '.-', color=PlotSpecDict['color'])
    
    # Add label to x axis.
    AxesHandle.set_xlabel('time relative to target entry (sec.)')
    
    # Add label to y axis.
    AxesHandle.set_ylabel(OrdinateName)
    
    # Title the axes object
    AxesHandle.set_title(OrdinateName + TitleLabel)
    
    # Add legend.
    plt.legend(loc='lower right')

def CumulativeWindowGen(BoundaryWindow, StepWidth, Direction):
    
    # This function generates a list of domains of increasing width over the 
    # entire, user-specified domain BoundaryWindow. Each cumulatively-growing
    # window in the list is one StepWidth longer than the one that precedes it
    # in the list.  The Direction arguement specifies whether windows 
    # cummulatively grow from the floor (positive) or ceiling (negative) of
    # the full domain.  
    
    # Generate a list of windows that grow positively from the floor of BoundaryWindow
    # if Direction is positive.
    if (Direction == 'positive'):
        
        # Generate the ceiling of each domain in list.
        ForwardArray = np.round(np.arange(BoundaryWindow[0] + StepWidth, 
                                  BoundaryWindow[1] + 0.5*StepWidth, 
                                  StepWidth)/StepWidth)*StepWidth
        
        # Generate the common floor of each domain in list, combine with ceilings
        # and transpose into an N by 2 array.
        ArrayOfWindows = np.array([BoundaryWindow[0]*np.ones_like(ForwardArray), 
                                      ForwardArray]).transpose()
    elif (Direction == 'negative'):
        
        # Generate the floor of each domain in list.
        BackwardArray = np.arange(BoundaryWindow[0], BoundaryWindow[1], 
                                  StepWidth)[::-1]
        
        # Generate the common ceiling of each domain in list, combine with floors
        # and transpose into an N by 2 array.
        ArrayOfWindows = np.array([BackwardArray, BoundaryWindow[1]*np.ones_like(
                                    BackwardArray)]).transpose()
    
    return ArrayOfWindows
                               
def SlidingWindowGen(BoundaryWindow, StepWidth, WindowWidth):
    
    # This simple function generates a list of arrays of domain boundaries
    # that are slid over the wide, user-specified domain: BoundaryWindow.
    # Windows are WindowWidth wide, and slid in increments of StepWidth.
    
    # Generate a list of domain lower bounds.
    LowerBounds = np.arange(BoundaryWindow[0], BoundaryWindow[1] - WindowWidth + StepWidth, 
                            StepWidth)
    
    # Generate a list of domain upper bounds.
    UpperBounds = np.arange(BoundaryWindow[0] + WindowWidth, BoundaryWindow[1] + StepWidth,
                            StepWidth)
    
    # Combine corresponding lower and upper bounds into a 2-D array and transpose.
    ArrayOfWindows = np.array([LowerBounds, UpperBounds]).transpose()
    
    return ArrayOfWindows

def PeriEventExtractorDictJoiner(PE_dict1,  PE_dict2):

    # This function joins together two PeriEventExtractorDict dictionaries 
    # produced by PeriEventExtractor_Trace fucnction (above) and called by
    # SlidingWindowAnalysisFunc.py. This functionality can be used
    # generate a combined dataset on which to run subsequent decoding analyses.
      
    (D1_NumTraces, D1_NumSamples) = PE_dict1['PEA_Array'].shape
    (D2_NumTraces, D2_NumSamples) = PE_dict2['PEA_Array'].shape
    
    # List of keys pointing to dictionary elements to join.
    ListOfToJoinKeys = ['PEA_Array', 'TargetsVec', 'TrialIndicesByEventDict']

    # Initialize the combined target dictionary
    CombinedDict = PE_dict1
    
    # Attempt to join peri-event activity arrays from the two dictionaries.  To
    # be successfully joined, the two arrays must have the same number of columns.
    # For that reason, joining operations are placed in error-exception blocks here.
    try:
        # Join PEA_Array elements
        CombinedDict['PEA_Array'] = np.vstack([PE_dict1['PEA_Array'], PE_dict2['PEA_Array']])

    except ValueError:
        
        print('Error: PEA_Array fields from the two dicts likely have different number of columns.')
        
    try:
        # Join TargetsVec elements
        CombinedDict['TargetsVec'] = np.vstack([PE_dict1['TargetsVec'], PE_dict2['TargetsVec']])

    except ValueError:
        
        print('Error: TargetVec fields from the two dicts likely have different number of columns.')
        
    # Join TrialIndicesByEvent sub-dicts.  Iterate through entires in RefEventsList
    # joining like entries from  both dicts as well as adjusting corresponding
    # referencing indices.
    for RefEvent in CombinedDict['RefEventsDict']['RefEventsList']:
        
        # Build index list by adding the number of elements in the first set
        # to each index of the second set, and then concatenatethe adjusted 
        # second set onto the first.
        IndicesToAppend = PE_dict2['TrialIndicesByEventDict'][RefEvent] + \
            D1_NumTraces*np.ones_like(PE_dict2['TrialIndicesByEventDict'][RefEvent])
        
        CombinedDict['TrialIndicesByEventDict'][RefEvent] = \
            np.hstack([CombinedDict['TrialIndicesByEventDict'][RefEvent], 
                      IndicesToAppend])
            
    return CombinedDict

def UniquePairsGenerator(CellIDs):
    
    # Make sure that the datatype of CellIDs list is a numpy array
    CellIDs = np.array(CellIDs, dtype=object)
    
    # Check that each entry in CellIDs is unique.  If not, tell user and then
    # remove duplicates.
    if (len(CellIDs) != len(np.unique(CellIDs))):
        
        print("Warning: input argument list contains duplicate entries.\nFiltering list before proceeding...")
        CellIDs = np.unique(CellIDs)
        
    # Count number of elements in list
    (nIDs,) = CellIDs.shape
    
    # Verify that 
    # Generate the sequential list of indices used for later iterations.
    IndicesList = np.arange(0, nIDs)
    
    # Initialize array to contain full list of pairs
    IndexPairsList = np.empty((nIDs*nIDs, 2), dtype=int)
    
    # Populate 2D array of all possible pairs.  Note that we are working with 
    # list entry indices here rather than the entries themselves.  This is 
    # because numpy's unique method can only sort an array-of-arrays if
    # array entries are numerical--it will not do nested sorts if array entries
    # are a non-numeric datatype.  Secondly, using indexes rather than labels
    # assumes that each entry in the CellIDs array is unique, which should be
    # the case for cell identity labels.  Recall uniqueness is checked above.
    for i in IndicesList:
        
        for j in IndicesList:
            
            IndexPairsList[i*nIDs+j] = np.array([i, j])
            
    ### Ignoring order of entries in each pair, keep only unique pairs.
    # First  sort entries in order of their last (second) dimension.
    IndexPairsList.sort(axis=-1)
    
    # For processing speed, write PairsList into a C-based ordering that is 
    # contiguous in memory.
    cIndexPairsList = np.ascontiguousarray(IndexPairsList)
    
    # Select unique entries    
    UniqueIndexPairsList, inv, ct = np.unique(cIndexPairsList, return_inverse=True,
                                         return_counts=True, axis=0)
    
    # Filter out subarrays with same value
    KeepFilt = (UniqueIndexPairsList[:,0] != UniqueIndexPairsList[:,1])
    UniqueIndexPairsList = UniqueIndexPairsList[KeepFilt, :]
    
    UniquePairsList = np.array([CellIDs[UniqueIndexPairsList[:, 0]], 
                                CellIDs[UniqueIndexPairsList[:, 1]]]).transpose()
    
    return UniqueIndexPairsList, UniquePairsList
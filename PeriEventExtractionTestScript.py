#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 22:07:25 2019

@author: thugwithyoyo
"""
import numpy as np
import CalciumImagingBehaviorAnalysis as CIBA
from collections import defaultdict
import CalciumImagingFluorProcessing as CIFP
import PartialLeastSquares as PLS
import hwfunclib as HFL


PathToBehavFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_both_001_1-01Dat.csv'
PathToFluorFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_001_Left_CellTraces.csv'

#RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts', 'M6CT_Exit_ts']
RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Set parameters for peri-event extraction
BoundaryWindow = [-1.0, 1.0]

# Set parameters for PLS
nLatents = 5
    
def PeriEventExtractor_Calcium(PathToBehavFile, PathToFluorFile, RefEventsList, BoundaryWindow):
    
    EventIndices = np.arange(0, len(RefEventsList))
    
    BehavDict = CIBA.CinePlexCSV_parser(PathToBehavFile)
        
    # f calcium transients and write into Pandas
    # dataframe.
    CellFluorTraces_Frame = CIFP.IDPS_CellFluorTracesParser(PathToFluorFile)
    
    # Initialize dict for storing peri-event activity; one key per reference event.
    PeriEventActivityArraysDict = defaultdict()
    
    # Initialize input argument arrays for Partial Least Squares
    PEA_Array = np.array([[]]).transpose()
    TargetsVec = np.array([[]], dtype=int).transpose()
    
    # Generate a dictionary of lists of indices for referencing trials for each
    # event.
    TrialIndicesByEventDict = defaultdict()
    nTrialsPerEvent = np.empty((len(RefEventsList),), dtype=int)
    
    # Populate PLS input arrays by stacking vertically sets of peri-event activity
    # snippets for each event type.
    for i in EventIndices:
        
        # Extract peri-event snippets about each event of the current event type. 
        PeriEventActivityArraysDict[RefEventsList[i]] = CIFP.PeriEventTraceExtractor(
                CellFluorTraces_Frame, BehavDict[RefEventsList[i]].values, 
                BoundaryWindow)
        
        # Determine the number of events in the current set.
        (nTrialsPerEvent[i], nSamplesPerTrace) = PeriEventActivityArraysDict[RefEventsList[i]].shape
    
        # If this iteration is the first, then initialize empty arrays of appropriate
        # size for stacking.
        if i == 0:
            
            PEA_Array = np.empty((0, nSamplesPerTrace))
            TargetsVec = np.empty((0, 1))
        
        # Aquire number of total trials to generate index list for current set.
        (nTotalTrials, nSamplesPerTrace) = PEA_Array.shape    
        
        # Generate and store list of indices for current event set that indicate
        # trace row locations in the composite activity array PEA_Array. 
        TrialIndicesByEventDict[RefEventsList[i]] = np.arange(nTotalTrials, 
                               nTotalTrials + nTrialsPerEvent[i])
        
        # Append current peri-activity set to the composite array.        
        PEA_Array = np.vstack((PEA_Array, PeriEventActivityArraysDict[RefEventsList[i]]))
        
        # Each event type will be assigned its own unique integer.  Stack sets of 
        # event integer lables in the same manner as the PEA_Array above.
        RefEventIntVec = i*np.ones((nTrialsPerEvent[i], 1), dtype=int)
        
        TargetsVec = np.vstack((TargetsVec, RefEventIntVec))
        
    OutputDict = defaultdict()
    OutputDict['PEA_Array'] = PEA_Array
    OutputDict['TargetsVec'] = TargetsVec
    OutputDict['RefEventsList'] = RefEventsList
    OutputDict['TrialIndicesByEventDict'] =  TrialIndicesByEventDict


def PLS_DecoderPerformance(PeriEventExtractorDict, nLatents, **kwargs):
    
    # Unpack dictionary contents
    PeriEventExtractorDict['PEA_Array'] = PEA_Array
    PeriEventExtractorDict['TargetsVec'] = TargetsVec
    PeriEventExtractorDict['RefEventsList'] = RefEventsList
    PeriEventExtractorDict['TrialIndicesByEventDict'] =  TrialIndicesByEventDict
    
    # Index number of included events
    EventIndices = np.arange(0, len(RefEventsList))
    
    # Initialize confusion matrix from which mutual information will be determined.
    ConfusionMatrix = np.zeros((len(RefEventsList), len(RefEventsList)), dtype=int)
        
    #### LEAVE ONE OUT DECODER PERFORMANCE EVALUATION
    # Iterate through each of the input traces in PEA_Array.  Select the current
    # trace to be the test trace and use the rest to compute the PLS model decoder
    # matrices (B, B_0).  Using the test trace make a prediction and record the 
    # outcome in the confusion matrix.
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    PEA_ArrayRowIndices = np.arange(0, nTotalTrials)
 
    
    if np.isin('trials_to_include', np.array(list(kwargs.keys()))):
        
        TrialSetIndices = kwargs['trials_to_include']
    
    else:
        
        TrialSetIndices = np.arange(0, nTotalTrials)
        
    # Iterate through each single trial in the peri-event activity array.  Find
    # the event type to which the current trial belongs and compute the regression
    # algorithm using the remaining (n - 1) set of trials.        
    for TestTrialIndex in TrialSetIndices:
    
        # Find the event pool that contains the current trial
        for Event in RefEventsList:
            
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
        B, B_0 = PLS.PLS1(PEA_Array[ComplementFilt,:], TargetsVec[ComplementFilt,:], 
                          nLatents)
        
    
        # Using the resulting PLS model, predict the event during which the test
        # trace was recorded.
        PredictedEventIndex = int(np.round(TestActivity @ B + B_0)[0,0])
            
        TrueEventIndex = EventIndices[RefEventsList == np.array(CurrentEvent)][0]
        
        ConfusionMatrix[TrueEventIndex, PredictedEventIndex] += 1
        
    # Compute the PLS regression matrices for the entire set of trials.    
    B, B_0 = PLS.PLS1(PEA_Array, TargetsVec, nLatents)
    
    #Predicted = X @ B + B_0
    Predicted = PEA_Array @ B + B_0
    
    return {    
            'confusion_matrix' : ConfusionMatrix,
            'performance' : HFL.PerformanceCalculator(ConfusionMatrix),
            'mutual_info' : HFL.MutualInformationFromConfusionMatrix(ConfusionMatrix),
            'targets' : TargetsVec,
            'predicted' : Predicted,
            'B' : B,
            'B_0' : B_0
            }
    
    #nIncorrects = np.sum(np.abs(np.squeeze(np.round(Predicted)) - TargetsVec.transpose()[0]))

##### SHUFFLE CONTROL
# To test valitidity of decoding, shuffle event lables prior to running decoder
#TargetsVec = TargetsVec.transpose()[0]
#np.random.shuffle(TargetsVec)
#TargetsVec = np.array([TargetsVec]).transpose()
    
##### PLS_MonteCarlo
def PLS_MonteCarlo(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel):
    
    # Calculate alphas and boundary indices from arguments
    alphas = np.empty((2,), dtype=int)
    alphas[0] = np.ceil(nRepetitions*(1. - ConfLevel)/2.)
    alphas[1] = np.floor(nRepetitions*(1. + ConfLevel)/2.)
    
    # Unpack needed dictionary contents
    PeriEventExtractorDict['PEA_Array'] = PEA_Array
    #PeriEventExtractorDict['TargetsVec'] = TargetsVec
    PeriEventExtractorDict['RefEventsList'] = RefEventsList
    #PeriEventExtractorDict['TrialIndicesByEventDict'] =  TrialIndicesByEventDict

    # Initialize arrays for storing simulated output distributions    
    nRefEvents = len(RefEventsList)
    ConfusionMatrixSimDist = np.empty((nRefEvents, nRefEvent, nRepetitions),
                                      dtype=int)
    PerformanceSimDist = np.empty((nRepetitions,), dtype=float)
    MutualInfoSimDist = np.empty((nRepetitions,), dtype=float)
    
    # Generate set of indices from which to extract sets for each bootstrap 
    # iteration
    (nTotalTrials, nTotalFeatures) = PEA_Array.shape
    #MasterIndicesSet = np.arange(0, nTotalTrials)
    
    for i in np.arange(0, nRepetitions):
        
        InclusionSet = np.random.randint(0, high=nTotalTrials, 
                                         size=(nTotalTrials,))
        
        PerformanceStruct = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    nLatents, trials_to_include=InclusionSet)
        
        ConfusionMatrixSimDist[:,:,i] = PerformanceStruct['confusion_matrix']
        PerformanceSimDist[i] = PerformanceStruct['performance']
        MutualInfoSimDist[i] = PerformanceStruct['mutual_info']
        
        
    return {
            'confusion_matrix_CLs' : np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas],
            'performance_CLs' : np.sort(PerformanceSimDist, axis=-1)[alphas],
            'mutual_info_CLs' : np.sort(MutualInfoSimDist, axis=-1)[alphas]
            }
    
PeriEventExtractorDict = PeriEventExtractor_Calcium(PathToBehavFile, PathToFluorFile, 
                           RefEventsList, BoundaryWindow)

PerformanceDict = PLS_DecoderPerformance(PeriEventExtractorDict, nLatents)
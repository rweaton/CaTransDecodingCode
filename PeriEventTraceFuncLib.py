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
import os
import shelve
from tkinter import *
from tkinter import filedialog
import matplotlib.pyplot as plt

#filename, file_extension = os.path.splitext('/path/to/somefile.ext')

# Generate dictionary of reference events
def BehavDictGen(PathToBehavFile):
    
    # Determine particular subroutine to parse behavioral events.
    _, file_extension = os.path.splitext(PathToBehavFile)
    
    if (file_extension == '.csv'):
        
        BehavDict = CIBA.CinePlexCSV_parser(PathToBehavFile)
        
    elif (file_extension == '.json'):
        
        BehavDict = CIBA.CinePlexFromMatJSON_parser(PathToBehavFile)
        
    else:
        
        print('Behavioral events: unrecognized file format.')
        
    return BehavDict

# Generate dataframe of calcium events
def CellFluorTraces_FrameGen(PathToFluorFile):
    
        # Determine particular subroutine to parse calcium traces.
    _, file_extension = os.path.splitext(PathToFluorFile)
    
    if (file_extension == '.csv'):
        
        CellFluorTraces_Frame = CIFP.IDPS_CellFluorTracesParser(PathToFluorFile)

    elif (file_extension == '.json'):
        
        CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)
        
    else:
        
        print('Calcium fluorescence traces: unrecognized file format.')

    return CellFluorTraces_Frame

def RemoveRepeatTargetEntries(BehavDict, RefEventsList, RelativeTolWindow):

    # Peripheral target entry events
    #RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

    #BehavDict = BehavDictGen(PathToBehavFile)
    nEventListsToFilter = len(RefEventsList)
    
    #RapidRepeatDict = np.empty((nEventListsToFilter,), dtype=dict)
    EventFilters = defaultdict(dict)
    
    for i in np.arange(0, nEventListsToFilter):
        
        #ComparatorDict = CIBA.EventComparator(tlist1, tlist2, tol_window)
        RapidRepeatDict =  CIBA.EventComparator(
                                        BehavDict[RefEventsList[i]].values, 
                                        BehavDict[RefEventsList[i]].values, 
                                        RelativeTolWindow)
        
        EventFilters[RefEventsList[i]] = CIBA.RemoveEventsInTolWindowFiltGen(
                                            RapidRepeatDict)

    return EventFilters
    
def PeriEventExtractor_Trace(BehavDict, CellFluorTraces_Frame, RefEventsDict, BoundaryWindow):
    
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
        B, B_0 = PLS.PLS1(PEA_Array[ComplementFilt,:], TargetsVec[ComplementFilt,:], 
                          nLatents)
        
    
        # Using the resulting PLS model, predict the event during which the test
        # trace was recorded.
        ModelOutput = TestActivity @ B + B_0
        
        EuclideanDiffVals = np.empty((EventIndices.size,), dtype=float)
        
        for i in EventIndices:
            
            EuclideanDiffVals[i] = HFL.EuclideanDistanceCalculator(np.array(
                    RefEventsDict['AssignedEventVals'][i]), ModelOutput)
        #EuclideanDiffVals = np.array(RefEventsDict['AssignedEventVals']) - ModelOutput
        
        PredictedEventIndex = np.argmin(EuclideanDiffVals)
        
#        PredictedEventIndex = int(np.round(TestActivity @ B + B_0)[0,0])
#        
#        if PredictedEventIndex > np.max(TargetsVec):
#            
#            PredictedEventIndex = np.max(TargetsVec)
#            
#        if (PredictedEventIndex < 0):
#            
#            PredictedEventIndex = 0
            
        TrueEventIndex = EventIndices[RefEventsDict['RefEventsList'] == np.array(CurrentEvent)][0]
        
        ConfusionMatrix[TrueEventIndex, PredictedEventIndex] += 1
        
    # Compute the PLS regression matrices for the entire set of trials.    
    B, B_0 = PLS.PLS1(PEA_Array, TargetsVec, nLatents)
    
    ModelOutput = PEA_Array @ B + B_0
    
    Predicted = np.empty((nTotalTrials,1), dtype=int)
    
    for j in np.arange(0, nTotalTrials):
        
        for i in EventIndices:
            
            EuclideanDiffVals[i] = HFL.EuclideanDistanceCalculator(np.array(
                    RefEventsDict['AssignedEventVals'][i]), ModelOutput[j])    
            
        PredictedEventIndex = np.argmin(EuclideanDiffVals)
            
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

    
##### PLS_MonteCarlo
def PLS_MonteCarlo(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel):
    
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
    
    for i in np.arange(0, nRepetitions):
        
        InclusionSet = np.random.randint(0, high=nTotalTrials, 
                                         size=(nTotalTrials,))
        
        PerformanceStruct = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                    nLatents, trials_to_include=InclusionSet)
        
        ConfusionMatrixSimDist[:,:,i] = PerformanceStruct['confusion_matrix']
        PerformanceSimDist[i] = PerformanceStruct['performance']
        MutualInfoSimDist[i] = PerformanceStruct['mutual_info']
        
        
    return {
            'confusion_matrix_median' : np.median(ConfusionMatrixSimDist, axis=-1),
            'confusion_matrix_CLs' : np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas],
            'performance_median' : np.median(PerformanceSimDist, axis=-1),
            'performance_CLs' : np.sort(PerformanceSimDist, axis=-1)[alphas],
            'mutual_info_median' : np.median(MutualInfoSimDist, axis=-1),
            'mutual_info_CLs' : np.sort(MutualInfoSimDist, axis=-1)[alphas]
            }
    
##### SHUFFLE CONTROL
# To test valitidity of decoding, shuffle event lables prior to running decoder
def PLS_Shuffle(PeriEventExtractorDict, nLatents, nRepetitions, ConfLevel):
    
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
        
    return {
            'confusion_matrix_median' : np.median(ConfusionMatrixSimDist, axis=-1),
            'confusion_matrix_CLs' : np.sort(ConfusionMatrixSimDist, axis=-1)[:,:,alphas],
            'performance_median' : np.median(PerformanceSimDist, axis=-1),
            'performance_CLs' : np.sort(PerformanceSimDist, axis=-1)[alphas],
            'mutual_info_median' : np.median(MutualInfoSimDist, axis=-1),
            'mutual_info_CLs' : np.sort(MutualInfoSimDist, axis=-1)[alphas]
            }
    
#TargetsVec = TargetsVec.transpose()[0]
#np.random.shuffle(TargetsVec)
#TargetsVec = np.array([TargetsVec]).transpose()
    
def ShelveWorkspace():
    
    ScriptDir = os.getcwd()
    root = Tk()
    root.filename =  filedialog.asksaveasfilename(
                        initialdir = "/home/thugwithyoyo/Desktop",
                        title = "Enter name of file to save", 
                        filetypes = (("out files","*.out"), ("all files","*.*")))    
    
    if root.filename is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    
    LoadFilePath = root.filename
    root.destroy()
    (PathToFile, Filename) = os.path.split(LoadFilePath)
    
    os.chdir(PathToFile)
     
    my_shelf = shelve.open(Filename)    
    #PathToSaveFile = root.filename
    #my_shelf = shelve.open(PathToSaveFile, 'n') # 'n' for new

    for key in dir():
        
        try:
            
            my_shelf[key] = globals()[key]
            
        except TypeError:
            #
            # __builtins__, my_shelf, and imported modules can not be shelved.
            #
            print('ERROR shelving: {0}'.format(key))
            
    #root.destroy()        
    my_shelf.close()
    
    os.chdir(ScriptDir)
    
def RestoreShelvedWorkspace():
    
    ScriptDir = os.getcwd()
    
    root = Tk()
    root.filename =  filedialog.askopenfilename(
            initialdir = "/home/thugwithyoyo/Desktop",
            title = "Select file to open",
            filetypes = (("dat files","*.dat"),("all files","*.*")))

    RestoreFilePath = root.filename
    root.destroy()
    (PathToFile, Filename) = os.path.split(RestoreFilePath)
    
    os.chdir(PathToFile)
     
    my_shelf = shelve.open(os.path.splitext(Filename)[0])
    
    for key in my_shelf:
        
        globals()[key]=my_shelf[key]
        
    my_shelf.close()
    
    os.chdir(ScriptDir)
    
def GenerateConfIntsPlot(ArrayOfConfIntsDicts, ArrayOfPerformanceDicts, PlotSpecDict, AxesHandle, Direction):
    
    (nDicts,) = ArrayOfConfIntsDicts.shape
    
    StepSize = np.max(np.abs(ArrayOfConfIntsDicts[1]['PeriEventDomain'] - 
                             ArrayOfConfIntsDicts[0]['PeriEventDomain']))
    
    X = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)
    Y_median = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)
    CI_Bars = np.empty((nDicts,2), dtype=float)
    CI_Bounds = np.empty((nDicts,2), dtype=float)
    Y_test = np.empty(ArrayOfConfIntsDicts.shape, dtype=float)
 
    for i in np.arange(0, nDicts):
        
        if (Direction == 'bw'):
            
            X[i] = np.min(ArrayOfConfIntsDicts[i]['PeriEventDomain']) + 0.5*StepSize
            TraceLabel = 'neg. cumulative span from '+ str(ArrayOfConfIntsDicts[i]['PeriEventDomain'][1]) +' sec.'
            
        elif (Direction == 'fw'):
            
            X[i] = np.max(ArrayOfConfIntsDicts[i]['PeriEventDomain']) - 0.5*StepSize
            TraceLabel = 'pos. cumulative span from '+ str(ArrayOfConfIntsDicts[i]['PeriEventDomain'][0]) +' sec.'
            
        elif (Direction == 'fw_sliding'):
            
            X[i] = np.mean(np.array(ArrayOfConfIntsDicts[i]['PeriEventDomain']))
            WindowWidth = np.diff(np.array(ArrayOfConfIntsDicts[0]['PeriEventDomain']))[0]
            TraceLabel = 'Window of width '+str(WindowWidth)+' sec. slid over ' \
                            +str(round(StepSize, 3))+' sec. increments.'
            
#        Y_median[i] = ArrayOfConfIntsDicts[i]['performance_median']
        Y_median[i] = ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']]
        
#        CI_Bars[i,:] = np.abs(ArrayOfConfIntsDicts[i]['performance_CLs'] - 
#                              ArrayOfConfIntsDicts[i]['performance_median'])
        CI_Bars[i,:] = np.abs(ArrayOfConfIntsDicts[i][PlotSpecDict['measure_CLs']] - 
                              ArrayOfConfIntsDicts[i][PlotSpecDict['measure_median']])        

#        CI_Bounds[i,:] = ArrayOfConfIntsDicts[i]['performance_CLs']
        CI_Bounds[i,:] = ArrayOfConfIntsDicts[i][PlotSpecDict['measure_CLs']]        
        
        
#        Y_test[i] = ArrayOfPerformanceDicts[i]['performance']
        Y_test[i] = ArrayOfPerformanceDicts[i][PlotSpecDict['measure']]
        
    #fig = plt.figure()
    
    #AxesHandle.errorbar(X, Y_median, yerr=CI_Bars.transpose(), lolims=True, uplims=True, label=TraceLabel)
    AxesHandle.fill_between(X, CI_Bounds.transpose()[0,:], CI_Bounds.transpose()[1,:], label=TraceLabel, alpha=0.7)
    AxesHandle.plot(X, Y_test, '.-')
    AxesHandle.set_xlabel('time relative to target entry (sec.)')
    AxesHandle.set_ylabel(PlotSpecDict['measure'])
    AxesHandle.set_title(PlotSpecDict['measure'] + ' dependence on span of activity window')
    
    plt.legend(loc='lower right')

def CumulativeWindowGen(BoundaryWindow, StepWidth, Direction):
    
    if (Direction == 'positive'):
        
        ForwardArray = np.round(np.arange(BoundaryWindow[0] + StepWidth, 
                                  BoundaryWindow[1] + 0.5*StepWidth, 
                                  StepWidth)/StepWidth)*StepWidth

        ArrayOfWindows = np.array([BoundaryWindow[0]*np.ones_like(ForwardArray), 
                                      ForwardArray]).transpose()
    elif (Direction == 'negative'):
        
        BackwardArray = np.arange(BoundaryWindow[0], BoundaryWindow[1], 
                                  StepWidth)[::-1]
        
        ArrayOfWindows = np.array([BackwardArray, BoundaryWindow[1]*np.ones_like(
                                    BackwardArray)]).transpose()
    
    return ArrayOfWindows
                               
def SlidingWindowGen(BoundaryWindow, StepWidth, WindowWidth):
    
    LowerBounds = np.arange(BoundaryWindow[0], BoundaryWindow[1] - WindowWidth + StepWidth, 
                            StepWidth)
    
    UpperBounds = np.arange(BoundaryWindow[0] + WindowWidth, BoundaryWindow[1] + StepWidth,
                            StepWidth)
    
    ArrayOfWindows = np.array([LowerBounds, UpperBounds]).transpose()
    
    return ArrayOfWindows
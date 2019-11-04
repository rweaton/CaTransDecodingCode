#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:31:49 2019

@author: thugwithyoyo
"""

import numpy as np
from PeriEventTraceFuncLib import *
from collections import defaultdict

# Paths to data in JSON formatted files
PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_unique_B.json'
PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_unique_C.json'
#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14-11-01-41_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14-11-01-41_C.json'

SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56_unique_AccumulatingWindow'

# Define ParamsDict, the dictionary that contains the parameters for
# PLS decoding, Bootstrapping, Shuffle control processes.
# Peripheral target entry events
ParamsDict = defaultdict(dict)
ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']


# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set parameters for peri-event extraction
ParamsDict['BoundaryWindow'] = [-1., 2.]
ParamsDict['StepWidth'] = 0.1
ParamsDict['WindowWidth'] = 0.4

# Set parameters for PLS
ParamsDict['NumLatents'] = 5

# Set parameters for Monte Carlo estimation of  confidence intervals
ParamsDict['NumRepetitions'] = 5
ParamsDict['ConfLevel'] = 0.95

# Specified anti-tolerance window, relative to target entry, for detecting and
# removing repeat entries that followed shortly after the initial entry.
ParamsDict['RelativeTolWindow'] = (0.0001, 2.5)
    
def VariableWindowAnalysisFunc(PathToBehavFile, PathToFluorFile, SavePath, ParamsDict):
    
    # Pack into a dict the two lists above that specify info for reference 
    # events
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
    
    ArrayOfWindows_fw = CumulativeWindowGen(ParamsDict['BoundaryWindow'], 
                                            ParamsDict['StepWidth'], 
                                            'positive')
    
    ArrayOfWindows_bw = CumulativeWindowGen(ParamsDict['BoundaryWindow'], 
                                            ParamsDict['StepWidth'], 
                                            'negative')
    
    # Generate the unfiltered behavior dictionary.
    BehavDict = BehavDictGen(PathToBehavFile)
    
    # Detect rapid repeats within each event list.
    EventFilters = RemoveRepeatTargetEntries(BehavDict, 
                                             RefEventsDict['RefEventsList'], 
                                             ParamsDict['RelativeTolWindow'])
    
    # Remove repeat events
    for ef in EventFilters:
        
        BehavDict[ef] = BehavDict[ef][EventFilters[ef]]
    
    # Generate the data frame of calcium transients.
    CellFluorTraces_Frame = CellFluorTraces_FrameGen(PathToFluorFile)
    
    # Grow window forwards from floor
    (NumDomains_fw, _) = ArrayOfWindows_fw.shape
    
    # Initialize an empty array to contain output dictionaries from the 
    # decoder cross-validation perfomance and monte carlo bootstrap routines
    Performance_fw = np.empty((NumDomains_fw,), dtype=dict)
    ConfInts_fw = np.empty((NumDomains_fw,), dtype=dict)
    EventsShuffled_fw = np.empty((NumDomains_fw,), dtype=dict)
    
    for i in np.arange(0, NumDomains_fw):
     
        PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                        CellFluorTraces_Frame, RefEventsDict, 
                                        ArrayOfWindows_fw[i])
        
        # Generate a set of indices to test the inclusion portion of the performance code.
        PEA_Array = PeriEventExtractorDict['PEA_Array']
        (NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
        InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))
    
        Performance_fw[i] = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                                   ParamsDict['NumLatents'])
        Performance_fw[i].update({'PeriEventDomain': ArrayOfWindows_fw[i]})
    
        ConfInts_fw[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        ConfInts_fw[i].update({'PeriEventDomain': ArrayOfWindows_fw[i]})
        
        EventsShuffled_fw[i] = PLS_Shuffle(PeriEventExtractorDict, 
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        EventsShuffled_fw[i].update({'PeriEventDomain': ArrayOfWindows_fw[i]})

        
    # Grow window backwards from ceiling
    (NumDomains_bw, _) = ArrayOfWindows_bw.shape
    
    # Initialize an empty array to contain output dictionaries from the 
    # decoder cross-validation routine.
    Performance_bw = np.empty((NumDomains_bw,), dtype=dict)
    ConfInts_bw = np.empty((NumDomains_bw,), dtype=dict)
    EventsShuffled_bw = np.empty((NumDomains_bw,), dtype=dict)
    
    for i in np.arange(0, NumDomains_bw):
     
        PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                        CellFluorTraces_Frame, RefEventsDict, 
                                        ArrayOfWindows_bw[i])
        
        PeriEventExtractorDict.update({'PeriEventDomain': ArrayOfWindows_bw[i]})
    
        # Generate a set of indices to test the inclusion portion of the performance code.
        PEA_Array = PeriEventExtractorDict['PEA_Array']
        (NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
        InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))
    
        Performance_bw[i] = PLS_DecoderPerformance(PeriEventExtractorDict, 
                                                   ParamsDict['NumLatents'])
        Performance_bw[i].update({'PeriEventDomain': ArrayOfWindows_bw[i]})
    
        ConfInts_bw[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        ConfInts_bw[i].update({'PeriEventDomain': ArrayOfWindows_bw[i]})
        
        EventsShuffled_bw[i] = PLS_Shuffle(PeriEventExtractorDict, 
                                        ParamsDict['NumLatents'], 
                                        ParamsDict['NumRepetitions'], 
                                        ParamsDict['ConfLevel'])
        EventsShuffled_bw[i].update({'PeriEventDomain': ArrayOfWindows_bw[i]})
        
    
    #####################################
    ########   Start shelving   #########
    #####################################
    
    # Assemble save path
    # get root directory of save path from path to calcium data
    # SavePath is an arguement above.  The following script requires it 
    exec(open('./ShelveWorkspaceScript.py').read())

VariableWindowAnalysisFunc(PathToBehavFile, PathToFluorFile, SavePath, ParamsDict)

RestoreFilePath = SavePath+'.dat'
exec(open('./RestoreShelvedWorkspaceScript.py').read())

#### Plot outcome measures #####
        
# Specify outcome measures to  be plotted
PerformancePlotSpecDict = {'measure': 'performance',
                       'measure_median': 'performance_median',
                       'measure_CLs': 'performance_CLs',
                       'color':'blue'}

ShuffledPerformancePlotSpecDict = {'measure': 'performance_median',
                       'measure_median': 'performance_median',
                       'measure_CLs': 'performance_CLs',
                       'color':'lightblue'}

MutInfoPlotSpecDict = {'measure': 'mutual_info',
                       'measure_median': 'mutual_info_median',
                       'measure_CLs': 'mutual_info_CLs',
                       'color':'blue'}

ShuffledMutInfoPlotSpecDict = {'measure': 'mutual_info_median',
                       'measure_median': 'mutual_info_median',
                       'measure_CLs': 'mutual_info_CLs',
                       'color':'lightblue'}
    
# Plot performance dependence on increasing peri-event window span
fig1, axs1 = plt.subplots()
#fig1.suptitle(PerformancePlotSpecDict['measure'])

# Plot forward-going accumulation performance 
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'blue'}

GenerateConfIntsPlot(ConfInts_fw, Performance_fw, PlotSpecDict, 
                     axs1, 'fw')

PlotSpecDict = {'measure': 'performance_median',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'lightblue'}

GenerateConfIntsPlot(EventsShuffled_fw, EventsShuffled_fw, PlotSpecDict, 
                     axs1, 'fw')

# Plot backward-going accumulation performance
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'green'}

GenerateConfIntsPlot(ConfInts_bw, Performance_bw, PlotSpecDict, 
                     axs1, 'bw')

PlotSpecDict = {'measure': 'performance_median',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'lightgreen'}

GenerateConfIntsPlot(EventsShuffled_bw, EventsShuffled_bw, PlotSpecDict, 
                     axs1, 'bw')

axs1.set_xbound(lower=ParamsDict['BoundaryWindow'][0], upper=ParamsDict['BoundaryWindow'][1])
axs1.set_ybound(lower=0.4, upper=1.)

# Plot mutual information dependence on increasing peri-event window span
fig2, axs2 = plt.subplots()
#fig2.suptitle(MutInfoPlotSpecDict['measure'])
PlotSpecDict = {'measure': 'mutual_info',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'color':'blue'}

GenerateConfIntsPlot(ConfInts_fw, Performance_fw, PlotSpecDict, 
                     axs2, 'fw')

PlotSpecDict = {'measure': 'mutual_info_median',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'color':'lightblue'}

GenerateConfIntsPlot(EventsShuffled_fw, EventsShuffled_fw, PlotSpecDict, 
                     axs2, 'fw')

PlotSpecDict = {'measure': 'mutual_info',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'color':'green'}

GenerateConfIntsPlot(ConfInts_bw, Performance_bw, PlotSpecDict, 
                     axs2, 'bw')

PlotSpecDict = {'measure': 'mutual_info_median',
                'measure_median': 'mutual_info_median',
                'measure_CLs': 'mutual_info_CLs',
                'color':'lightgreen'}

GenerateConfIntsPlot(EventsShuffled_bw, EventsShuffled_bw, PlotSpecDict, 
                     axs2, 'bw')

axs2.set_xbound(lower=ParamsDict['BoundaryWindow'][0], upper=ParamsDict['BoundaryWindow'][1])
axs2.set_ybound(lower=0., upper=1.)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:31:49 2019

@author: thugwithyoyo
"""

import numpy as np
from PeriEventTraceFuncLib import *

# Paths to data in JSON formatted files
#PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_B.json'
#PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_C.json'
PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14-11-01-41_B.json'
PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-14-11-01-41_C.json'

# Peripheral target entry events
RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
AssignedEventVals = [-1, 1]

# Pack into a dict the two lists above that specify info for reference 
# events
RefEventsDict = {'RefEventsList' : RefEventsList,
                 'AssignedEventVals' : AssignedEventVals}

# Specify outcome measures to  be plotted
PerformancePlotSpecDict = {'measure': 'performance',
                       'measure_median': 'performance_median',
                       'measure_CLs': 'performance_CLs'}

MutInfoPlotSpecDict = {'measure': 'mutual_info',
                       'measure_median': 'mutual_info_median',
                       'measure_CLs': 'mutual_info_CLs'}

# Set parameters for peri-event extraction
BoundaryWindow = [-1., 2.]
StepWidth = 0.1
WindowWidth = 0.4


#ForwardArray = np.round(np.arange(BoundaryWindow[0] + StepWidth, 
#                                  BoundaryWindow[1] + 0.5*StepWidth, 
#                                  StepWidth)/StepWidth)*StepWidth
#                                  
#BackwardArray = np.arange(BoundaryWindow[0], BoundaryWindow[1], 
#                          StepWidth)[::-1]
#
#ArrayOfWindows_fw = np.array([BoundaryWindow[0]*np.ones_like(ForwardArray), 
#                              ForwardArray]).transpose()
#
#ArrayOfWindows_bw = np.array([BackwardArray, BoundaryWindow[1]*np.ones_like(BackwardArray)]).transpose()

ArrayOfWindows_fw = CumulativeWindowGen(BoundaryWindow, StepWidth, 'positive')

ArrayOfWindows_bw = CumulativeWindowGen(BoundaryWindow, StepWidth, 'negative')

# Set parameters for PLS
NumLatents = 5

# Set parameters for Monte Carlo estimation of  confidence intervals
NumRepetitions = 30
ConfLevel = 0.95

# Specified anti-tolerance window, relative to target entry, for detecting and
# removing repeat entries that followed shortly after the initial entry.
RelativeTolWindow = (0.0001, 2.5)

# Generate the unfiltered behavior dictionary.
BehavDict = BehavDictGen(PathToBehavFile)

# Detect rapid repeats within each event list.
EventFilters = RemoveRepeatTargetEntries(BehavDict, RefEventsList, 
                                         RelativeTolWindow)

# Remove repeat events
for ef in EventFilters:
    
    BehavDict[ef] = BehavDict[ef][EventFilters[ef]]

CellFluorTraces_Frame = CellFluorTraces_FrameGen(PathToFluorFile)

# Grow window forwards from floor
(NumDomains_fw, _) = ArrayOfWindows_fw.shape

# Initialize an empty array to contain output dictionaries from the 
# decoder cross-validation perfomance and monte carlo bootstrap routines
Performance_fw = np.empty((NumDomains_fw,), dtype=dict)
ConfInts_fw = np.empty((NumDomains_fw), dtype=dict)

for i in np.arange(0, NumDomains_fw):
 
    PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                    CellFluorTraces_Frame, RefEventsDict, 
                                    ArrayOfWindows_fw[i])
    
    # Generate a set of indices to test the inclusion portion of the performance code.
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    (NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
    InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))

    Performance_fw[i] = PLS_DecoderPerformance(PeriEventExtractorDict, NumLatents)
    Performance_fw[i].update({'PeriEventDomain': ArrayOfWindows_fw[i]})

    ConfInts_fw[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                               NumLatents, NumRepetitions, 
                                               ConfLevel)
    ConfInts_fw[i].update({'PeriEventDomain': ArrayOfWindows_fw[i]})
    
# Grow window backwards from ceiling
(NumDomains_bw, _) = ArrayOfWindows_bw.shape

# Initialize an empty array to contain output dictionaries from the 
# decoder cross-validation routine.
Performance_bw = np.empty((NumDomains_bw,), dtype=dict)
ConfInts_bw = np.empty((NumDomains_bw), dtype=dict)

for i in np.arange(0, NumDomains_bw):
 
    PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                    CellFluorTraces_Frame, RefEventsDict, 
                                    ArrayOfWindows_bw[i])
    
    PeriEventExtractorDict.update({'PeriEventDomain': ArrayOfWindows_bw[i]})

    # Generate a set of indices to test the inclusion portion of the performance code.
    PEA_Array = PeriEventExtractorDict['PEA_Array']
    (NumTotalTrials, NumTotalFeatures) = PEA_Array.shape
    InclusionSet = np.random.randint(0, high=NumTotalTrials, size=(NumTotalTrials,))

    Performance_bw[i] = PLS_DecoderPerformance(PeriEventExtractorDict, NumLatents)
    Performance_bw[i].update({'PeriEventDomain': ArrayOfWindows_bw[i]})

    ConfInts_bw[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                               NumLatents, NumRepetitions, 
                                               ConfLevel)
    ConfInts_bw[i].update({'PeriEventDomain': ArrayOfWindows_bw[i]}) 
    
#### Plot outcome measures #####
    
# Plot performance dependence on increasing peri-event window span
fig1, axs1 = plt.subplots()
#fig1.suptitle(PerformancePlotSpecDict['measure'])
GenerateConfIntsPlot(ConfInts_fw, Performance_fw, PerformancePlotSpecDict, 
                     axs1, 'fw')
GenerateConfIntsPlot(ConfInts_bw, Performance_bw, PerformancePlotSpecDict, 
                     axs1, 'bw')
axs1.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs1.set_ybound(lower=0.4, upper=1.)

# Plot mutual information dependence on increasing peri-event window span
fig2, axs2 = plt.subplots()
#fig2.suptitle(MutInfoPlotSpecDict['measure'])
GenerateConfIntsPlot(ConfInts_fw, Performance_fw, MutInfoPlotSpecDict, 
                     axs2, 'fw')
GenerateConfIntsPlot(ConfInts_bw, Performance_bw, MutInfoPlotSpecDict, 
                     axs2, 'bw')
axs2.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs2.set_ybound(lower=0., upper=1.)


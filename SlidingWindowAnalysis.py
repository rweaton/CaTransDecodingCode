#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 09:33:08 2019

@author: thugwithyoyo
"""
import numpy as np
from PeriEventTraceFuncLib import *

# Paths to data in JSON formatted files
PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_unique_B.json'
PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_unique_C.json'

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

ShuffledPerformancePlotSpecDict = {'measure': 'performance_median',
                       'measure_median': 'performance_median',
                       'measure_CLs': 'performance_CLs'}

MutInfoPlotSpecDict = {'measure': 'mutual_info',
                       'measure_median': 'mutual_info_median',
                       'measure_CLs': 'mutual_info_CLs'}

ShuffledMutInfoPlotSpecDict = {'measure': 'mutual_info_median',
                       'measure_median': 'mutual_info_median',
                       'measure_CLs': 'mutual_info_CLs'}

# Set parameters for peri-event extraction
BoundaryWindow = [-1., 2.]
StepWidth = 0.1
WindowWidth = 0.4

ArrayOfSlidWindows = SlidingWindowGen(BoundaryWindow, StepWidth, WindowWidth)

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
#EventFilters = RemoveRepeatTargetEntries(BehavDict, RefEventsList, 
#                                         RelativeTolWindow)

# Remove repeat events
#for ef in EventFilters:
    
#    BehavDict[ef] = BehavDict[ef][EventFilters[ef]]

# Generate the data frame of calcium transients.
CellFluorTraces_Frame = CellFluorTraces_FrameGen(PathToFluorFile)

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

    Performance[i] = PLS_DecoderPerformance(PeriEventExtractorDict, NumLatents)
    Performance[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})

    ConfInts[i] = PLS_MonteCarlo(PeriEventExtractorDict,
                                               NumLatents, NumRepetitions, 
                                               ConfLevel)
    ConfInts[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
    
    EventsShuffled[i] = PLS_Shuffle(PeriEventExtractorDict, 
                                              NumLatents, NumRepetitions, 
                                              ConfLevel)
    EventsShuffled[i].update({'PeriEventDomain': ArrayOfSlidWindows[i]})
    
    
#### Plot outcome measures #####
    
# Plot performance dependence on increasing peri-event window span
fig1, axs1 = plt.subplots()
#fig1.suptitle(PerformancePlotSpecDict['measure'])
GenerateConfIntsPlot(ConfInts, Performance, PerformancePlotSpecDict, 
                     axs1, 'fw_sliding')

GenerateConfIntsPlot(EventsShuffled, EventsShuffled, ShuffledPerformancePlotSpecDict, 
                     axs1, 'fw_sliding')
axs1.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs1.set_ybound(lower=0.4, upper=1.)

# Plot mutual information dependence on increasing peri-event window span
fig2, axs2 = plt.subplots()
#fig2.suptitle(MutInfoPlotSpecDict['measure'])
GenerateConfIntsPlot(ConfInts, Performance, MutInfoPlotSpecDict, 
                     axs2, 'fw_sliding')
GenerateConfIntsPlot(EventsShuffled, EventsShuffled, ShuffledMutInfoPlotSpecDict, 
                     axs2, 'fw_sliding')
axs2.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs2.set_ybound(lower=0., upper=1.)
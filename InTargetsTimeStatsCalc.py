#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:38:04 2019

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
#RefEventsList = ['M6T0_Entry_ts', 'M6T0_Exit_ts', 'M6T1_Entry_ts', 'M6T1_Exit_ts']
RefEventsList = ['M6T0_Exit_ts', 'M6T1_Exit_ts']

BehavDict = BehavDictGen(PathToBehavFile)

#RelativeTolWindow = (0.0001, 5.)
RelativeTolWindow = (-2.5, -0.0001)

# Remove rapid repeats
EventFilters = RemoveRepeatTargetEntries(BehavDict, RefEventsList, RelativeTolWindow)

BehavDict['M6T0_Exit_ts'] = BehavDict['M6T0_Exit_ts'][EventFilters['M6T0_Exit_ts']]
BehavDict['M6T1_Exit_ts'] = BehavDict['M6T1_Exit_ts'][EventFilters['M6T1_Exit_ts']]


RelativeTolWindow = (-5., -0.0001)
T0_InTimesDict =  CIBA.EventComparator(BehavDict['M6T0_Exit_ts'].values, 
                                       BehavDict['M6T0_Entry_ts'].values, 
                                       RelativeTolWindow)

T1_InTimesDict =  CIBA.EventComparator(BehavDict['M6T1_Exit_ts'].values, 
                                       BehavDict['M6T1_Entry_ts'].values, 
                                       RelativeTolWindow)

T0_InTargDurations = np.empty((T0_InTimesDict.size,))
T1_InTargDurations = np.empty((T1_InTimesDict.size,))

for i in np.arange(0, T0_InTimesDict.size):
    
    T0_InTargDurations[i] = BehavDict['M6T0_Exit_ts'].values[i] - \
                            T0_InTimesDict[i]['within_tol_ts'][-1]
                                        

for i in np.arange(0, T1_InTimesDict.size):
    
    T1_InTargDurations[i] = BehavDict['M6T1_Exit_ts'].values[i] - \
                            T1_InTimesDict[i]['within_tol_ts'][-1]
                                        
 
fig, axs = plt.subplots(nrows=1, ncols=2)
axs[0].hist(T0_InTargDurations, range=(0., 4.), bins=30)
axs[1].hist(T1_InTargDurations, range=(0., 4.), bins=30)    
axs[0].set_title('Right Target (T0)')
axs[1].set_title('Left Target (T1)') 
axs[0].set_xlabel('in-target duration (sec.)')
axs[1].set_xlabel('in-target duration (sec.)')
axs[0].set_ylim(0., max(axs[0].get_ylim()[-1], axs[1].get_ylim()[-1]))
axs[1].set_ylim(0., max(axs[0].get_ylim()[-1], axs[1].get_ylim()[-1]))

path = os.path.normpath(PathToBehavFile)
fig.suptitle('In-target durations ' + path.split(os.sep)[-1])

            
#EventFilter = CIBA.KeepOnlyEventsInTolWindowFiltGen(T0_InTimesDict)

#EventFilters = RemoveRepeatTargetEntries(BehavDict, RefEventsList, RelativeTolWindow)

#M6T0_IntertrialDeltaTs = np.diff(BehavDict['M6T0_Entry_ts'][T6T0_EventFilter])

#ComparatorDict = CIBA.EventComparator(tlist1, tlist2, tol_window)
#M6T0_RapidRepeatDict =  CIBA.EventComparator(BehavDict['M6T0_Entry_ts'].values, 
#                                             BehavDict['M6T0_Entry_ts'].values, 
#                                             (0.0001, 2.5))
#
#T6T0_EventFilter = CIBA.RemoveEventsInTolWindowFiltGen(M6T0_RapidRepeatDict)
#M6T0_IntertrialDeltaTs = np.diff(BehavDict['M6T0_Entry_ts'][T6T0_EventFilter])

#M6T1_RapidRepeatDict =  CIBA.EventComparator(BehavDict['M6T1_Entry_ts'].values, 
#                                             BehavDict['M6T1_Entry_ts'].values, 
#                                             (0.0001, 2.5))
#
#T6T1_EventFilter = CIBA.RemoveEventsInTolWindowFiltGen(M6T1_RapidRepeatDict)
#M6T1_IntertrialDeltaTs = np.diff(BehavDict['M6T1_Entry_ts'][T6T1_EventFilter])

#fig, axs = plt.subplots(nrows=1, ncols=2)
#axs[0].hist(M6T0_IntertrialDeltaTs, range=(0., 20.))
#axs[1].hist(M6T1_IntertrialDeltaTs, range=(0., 20.))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 15 23:55:00 2019

@author: thugwithyoyo
"""

import numpy as np
from PeriEventTraceFuncLib import *
import subprocess

def ShelvingTestor(SavePath):
    
    PathToBehavFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-10-11-37-56_unique_B.json'
    
    # Some  test variables to save.
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
    
    
    #####################################
    ########   Start shelving   #########
    #####################################
    
    #ShelveWorkspace(SavePath)
    exec(open('./ShelveWorkspaceScript.py').read())
    #subprocess.call('./ShelveWorkspaceScript.py')
    #import ShelveWorkspaceScript

SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/TestSave2'
ShelvingTestor(SavePath)
#RestorePath = SavePath+'.dat'
#RestoreShelvedWorkspace(RestorePath)
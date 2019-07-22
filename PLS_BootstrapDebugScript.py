#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 14:21:44 2019

@author: thugwithyoyo
"""
import numpy as np
from PeriEventTraceFuncLib import *

# Paths to data in .csv files
#PathToBehavFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_both_001_1-01Dat.csv'
#PathToFluorFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_001_Left_CellTraces.csv'

# Paths to data in JSON formatted files
#PathToBehavFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/2018-12-10-11-37-56_B.json'
#PathToFluorFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/2018-12-10-11-37-56_C.json'
PathToBehavFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/2018-12-05-11-29-45_B.json'
PathToFluorFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/2018-12-05-11-29-45_C.json'

#RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts', 'M6CT_Exit_ts']
#RefEventsList = ['M6T0_Entry_ts', 'M6T1_Entry_ts']
RefEventsList = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
AssignedEventVals = [-1, 1]
RefEventsDict = {'RefEventsList' : RefEventsList,
                 'AssignedEventVals' : AssignedEventVals}

# Set parameters for peri-event extraction
BoundaryWindow = [-1.0, 1.0]

# Set parameters for PLS
nLatents = 5

# Set parameters for Monte Carlo estimation of  confidence intervals
nRepetitions = 10
ConfLevel = 0.95

BehavDict = BehavDictGen(PathToBehavFile)

CellFluorTraces_Frame = CellFluorTraces_FrameGen(PathToFluorFile)

PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, CellFluorTraces_Frame, 
                                                  RefEventsDict, BoundaryWindow)

# Generate a set of indices to test the inclusion portion of the performance code.
PEA_Array = PeriEventExtractorDict['PEA_Array']
(nTotalTrials, nTotalFeatures) = PEA_Array.shape
InclusionSet = np.random.randint(0, high=nTotalTrials, size=(nTotalTrials,))

#PerformanceDict = PLS_DecoderPerformance(PeriEventExtractorDict, nLatents,
#                                               trials_to_include=InclusionSet)
# Run performance statistics on sample set.
PerformanceDict = PLS_DecoderPerformance(PeriEventExtractorDict, nLatents)

# Bootstrap to obtain medians and confidence intervals of performance measures
MonteCarlo_ConfIntsDict = PLS_MonteCarlo(PeriEventExtractorDict,
                                               nLatents, nRepetitions, 
                                               ConfLevel)

# Shuffle lables to evaluate chance band of decoder.
Shuffled_ConfIntsDict = PLS_Shuffle(PeriEventExtractorDict,
                                    nLatents, nRepetitions,
                                    ConfLevel)

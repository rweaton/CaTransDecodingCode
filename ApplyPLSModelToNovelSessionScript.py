#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 16:52:07 2020

@author: thugwithyoyo
"""
import numpy as np
import PeriEventTraceFuncLib as PETFL
import LongRegFuncLib as LRFL 

from collections import defaultdict

import os
import shelve
from tkinter import *
from tkinter import filedialog


# Specify the number of latent factors to include in PLS regression model.
nLatents = 3

############## Acquire the training and testing data sets #####################
# Acquire path of workspace to load.

# 
CellFluorTraces_Frame_train, BehavDict_train, SessionName_train = LRFL.BuildSessionDatasets()

#
CellFluorTraces_Frame_test, BehavDict_test, SessionName_test = LRFL.BuildSessionDatasets()

## Set defualt parent directories
#DataRootDir = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData'
#SaveRootDir = '/home/thugwithyoyo/CaTransDecoding/Output'
#
## Start file dialogs 
#root = tk.Tk()
#
### Begin Training block ##
## Prompt user to navigate to, and select, the session behavior file
#PathToBehavFile = askopenfilename(title='Select TRAINING session behavior file',
#                                  filetypes=[("json files","*.json"), 
#                                             ("csv files", "*.csv")],
#                                  initialdir=DataRootDir)
#
## Prompt user to navigate to, and select, the session imaging file
#PathToFluorFile = askopenfilename(title='Select TRAINING session fluorescence record file',
#                                  filetypes=[("json files","*.json"),
#                                             ("csv files", "*.csv")],
#                                  initialdir=DataRootDir)
#
## Assemble the dataframe of cell fluorescence traces.
#CellFluorTraces_Frame_train = CellFluorTraces_FrameGen(PathToFluorFile)
#
## Assemble the dictionary of event occurence timestamps
#BehavDict_train = BehavDictGen(PathToBehavFile)
### End Training block ##
#
### Begin Testing block ##
## Prompt user to navigate to, and select, the session behavior file
#PathToBehavFile = askopenfilename(title='Select TESTING session behavior file',
#                                  filetypes=[("json files","*.json"), 
#                                             ("csv files", "*.csv")],
#                                  initialdir=DataRootDir)
#
## Prompt user to navigate to, and select, the session imaging file
#PathToFluorFile = askopenfilename(title='Select TESTING session fluorescence record file',
#                                  filetypes=[("json files","*.json"),
#                                             ("csv files", "*.csv")],
#                                  initialdir=DataRootDir)
#
## Assemble the dataframe of cell fluorescence traces.
#CellFluorTraces_Frame_test = CellFluorTraces_FrameGen(PathToFluorFile)
#
## Assemble the dictionary of event occurence timestamps
#BehavDict_test = BehavDictGen(PathToBehavFile)
### End Testing block ##
#
## End file dialogs
#root.withdraw()


######### Arrange the fluorescence dataframe into global arrangement ##########

lr_frame = LRFL.BuildLongRegDataframeFromCSV()

#SessionName_test = '2018-12-27-11-04-59'
CommonCellsDict = LRFL.CommonCellsExtractor(lr_frame, 
                                            SessionName_train, 
                                            SessionName_test)

(NumCells,) = CommonCellsDict['S1_local_cellset'].shape

CellFluorTraces_Frame_train = CellFluorTraces_Frame_train[
        np.hstack([np.array(['Timestamps']), 
                   CommonCellsDict['S1_local_cellset']])]

CellFluorTraces_Frame_train_shuf = CellFluorTraces_Frame_train[
        np.hstack([np.array(['Timestamps']), 
                   np.random.choice(CommonCellsDict['S1_local_cellset'], 
                                    size=NumCells, replace=False)])]

CellFluorTraces_Frame_test = CellFluorTraces_Frame_test[
        np.hstack([np.array(['Timestamps']), 
                   CommonCellsDict['S2_local_cellset']])]

CellFluorTraces_Frame_test_shuf = CellFluorTraces_Frame_test[
        np.hstack([np.array(['Timestamps']), 
                   np.random.choice(CommonCellsDict['S2_local_cellset'], 
                                    size=NumCells, replace=False)])]

#"ColumnsSortMap_train" = pd.
# CellFlourTraces_Frame_train = CellFlourTraces_Frame_train["ColumnsSortMap_train"]

#"ColumnsSortMap_test" = pd.
# CellFlourTraces_Frame_test = CellFlourTraces_Frame_test["ColumnsSortMap_test"]

###############################################################################
# Run the peri event extractor using the (arranged) fluorescence dataframe and
# corresponding behavior dictionary to assemble the training and testing datasets.

# Define a dictionary to contain keyed access to sliding window analysis parameters
ParamsDict = defaultdict(dict)

# Right arm to both right and left targets
ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set peri-event extraction window (seconds)
ParamsDict['BoundaryWindow'] = [-1., 1.]

# Pack into a dict the two lists that specify the names and assigned 
# numerical values of the behavioral event types listed in ParamsDict
RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                 'AssignedEventVals' : ParamsDict['AssignedEventVals']}

# Run peri-event extraction and assemble the dictionary to be used for
# PLS-regression training.
PeriEventExtractorDict_train = PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                        CellFluorTraces_Frame_train, 
                                        RefEventsDict, 
                                        ParamsDict['BoundaryWindow'])

PeriEventExtractorDict_train_shuf = PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                        CellFluorTraces_Frame_train_shuf, 
                                        RefEventsDict, 
                                        ParamsDict['BoundaryWindow'])

# Run peri-event extraction and assemble the dictionary to be used to test the
# PLS-regression model above.
PeriEventExtractorDict_test = PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                        CellFluorTraces_Frame_test, 
                                        RefEventsDict, 
                                        ParamsDict['BoundaryWindow'])

PeriEventExtractorDict_test_shuf = PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                        CellFluorTraces_Frame_test_shuf, 
                                        RefEventsDict, 
                                        ParamsDict['BoundaryWindow'])

###############################################################################
# Compute the PLS regression matrices from the training dataset.

Performance_train = PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                           nLatents)

###############################################################################
## 

Performance_test = PETFL.ApplyPLSModel(Performance_train['B'], 
                                 Performance_train['B_0'], 
                                 PeriEventExtractorDict_test)

Performance_test_shuf = PETFL.ApplyPLSModel(Performance_train['B'], 
                                      Performance_train['B_0'], 
                                      PeriEventExtractorDict_test_shuf)

Performance_train_shuf = PETFL.ApplyPLSModel(Performance_train['B'], 
                                       Performance_train['B_0'], 
                                       PeriEventExtractorDict_train_shuf)
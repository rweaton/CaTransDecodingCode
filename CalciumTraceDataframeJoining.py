#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 04:31:24 2019

@author: thugwithyoyo
"""
import pandas as pd
import CalciumImagingFluorProcessing as CIFP
import numpy as np
from collections import defaultdict

PathToFile_Left='/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-27/2018-12-27-11-35-13_new_unique_C.json'
PathToFile_Right='/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-12-27/2018-12-27-11-35-22_new_unique_C.json'

FluorDataframe_Left = CIFP.FrameFromJSON(PathToFile_Left)
FluorDataframe_Right = CIFP.FrameFromJSON(PathToFile_Right)

def JoinDataFramesAcrossColumns(FluorDataframe_Left, FluorDataframe_Right):
    
    IndexLength = min([FluorDataframe_Left.shape[0], FluorDataframe_Right.shape[0]])
    
    nCols_left = FluorDataframe_Left.shape[1]
    nCols_right = FluorDataframe_Right.shape[1]
    
    Timestamps_df = pd.DataFrame(FluorDataframe_Left._slice(slice(0, 1), 1), 
                                 index=np.arange(0, IndexLength))
    
    FluorDataframe_Left = pd.DataFrame(FluorDataframe_Left._slice(slice(1, nCols_left), 1),
                                       index=np.arange(0, IndexLength))
    
    # Append "L" to column names here.
    NewColNames = []
    for ColName in list(FluorDataframe_Left):
        NewColNames.append(ColName+'L')
        
    FluorDataframe_Left.columns = NewColNames
                
    
    FluorDataframe_Right = pd.DataFrame(FluorDataframe_Right._slice(slice(1, nCols_right), 1), 
                                        index=np.arange(0, IndexLength))
    
    # Append "R" to column names here
    NewColNames = []
    for ColName in list(FluorDataframe_Right):
        NewColNames.append(ColName + 'R')
        
    FluorDataframe_Right.columns = NewColNames
        
    FluorDataframe_Combined = pd.concat([Timestamps_df, FluorDataframe_Left, FluorDataframe_Right], 
                                         axis=1, sort=False)
    
    return FluorDataframe_Combined

def RecordingJoiner(CaImag_df1, CaImag_df2, BehavDict1, BehavDict2):
    
    CaImag_df1_ts = pd.DataFrame(CaImag_df1._slice(slice(0, 1), 1)).values
    CaImag_df2_ts = pd.DataFrame(CaImag_df2._slice(slice(0, 1), 1)).values
    
    CaImag_df2_ts_offset = CaImag_df1_ts[-1, :]

    CombinedCaImag_df_ts = np.vstack([CaImag_df1_ts, 
                                CaImag_df2_ts + CaImag_df2_ts_offset])
    
    CombinedCaImag_df = pd.concat([CaImag_df1, CaImag_df2], axis=0)
    CombinedCaImag_df['Timestamps'] = CombinedCaImag_df_ts
            
    Filt = np.char.endswith(np.array(list(BehavDict2.keys())),'ts')
    ts_KeysList = np.array(list(BehavDict2.keys()))[Filt]
        
    CombinedBehavDict = defaultdict(dict)
    
    for Event in ts_KeysList:
        
        CombinedBehavDict[Event] = \
            pd.Series(data=np.hstack([BehavDict1[Event].values, 
                    BehavDict2[Event].values + CaImag_df2_ts_offset]))
    
    return CombinedCaImag_df, CombinedBehavDict
    
def TraceSampler(CellFluorTraces_Frame, nTracesToSample):
    
    # Acquire dataframe column names. 
    ColNames = np.array(CellFluorTraces_Frame.columns)

    # Randomly sample columns from date frame.    
    # Determine total number of columns in dataframe.  The first column 
    # contains corresponding timestamps.
    (IndexLength, nCols) = CellFluorTraces_Frame.shape
    
    # Create a list of integers that index the trace columns.
    TracesIndices = np.arange(1, nCols)

    # Randomly sample, without replacement, from the above list.  Append zero
    # to list so that timestamps column will be extracted as well.     
    SamplingSetIndices = np.hstack([np.array([0]), np.random.choice(TracesIndices, 
                                    nTracesToSample, replace=False)])
    
    # Return extracted dataframe subset.
    return  CellFluorTraces_Frame[ColNames[SamplingSetIndices]]

FluorDataframe_Combined = JoinDataFramesAcrossColumns(FluorDataframe_Left, FluorDataframe_Right)
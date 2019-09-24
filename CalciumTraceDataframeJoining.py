#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 04:31:24 2019

@author: thugwithyoyo
"""
import pandas as pd
import CalciumImagingFluorProcessing as CIFP
import numpy as np

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

FluorDataframe_Combined = JoinDataFramesAcrossColumns(FluorDataframe_Left, FluorDataframe_Right)
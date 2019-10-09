#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  6 15:49:00 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd
from PeriEventTraceFuncLib import *
from CalciumTraceDataframeFuncLib import *
from collections import defaultdict


PathsToBehavFiles = ['/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_B.json',
                     '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_B.json']

PathsToFluorFiles = ['/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-36-02_new_unique_C.json',
                     '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2019-01-24/2019-01-24-11-50-23_new_unique_C.json']


CaImag_df_Combined = pd.DataFrame([])
BehavDict_Combined = defaultdict(dict)

nFiles = len(PathsToBehavFiles)

for i in np.arange(0,nFiles):
    
    if i == 0:
        
        CaImag_df_Combined = CellFluorTraces_FrameGen(PathsToFluorFiles[i])
        BehavDict_Combined = BehavDictGen(PathsToBehavFiles[i])
        
    else:
 
        CaImag_df = CellFluorTraces_FrameGen(PathsToFluorFiles[i])
        BehavDict = BehavDictGen(PathsToBehavFiles[i])
        
        CaImag_df_Combined, BehavDict_Combined = \
            RecordingJoiner(CaImag_df_Combined, CaImag_df, 
                            BehavDict_Combined, BehavDict)
            
    
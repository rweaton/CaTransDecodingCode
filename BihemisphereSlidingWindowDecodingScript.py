#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 08:13:05 2020

@author: thugwithyoyo
"""

import CalciumTraceDataframeFuncLib as CTDFL
import SlidingWindowAnalysisFunc as SWAF

import os
import shelve

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.messagebox import askyesno

#RestoreFilePath = SavePath +'.dat'

# Have user select workspace to load.  Workspace must be a shelve object '.dat'
# file generated from DecodingStabilityBatchProcessorFunc.py
DefaultWorkspaceDir = '/home/thugwithyoyo/CaTransDecoding/Output'

root = tk.Tk()

RestoreFilePath = askopenfilename(initialdir=DefaultWorkspaceDir, 
                                  title="Select a shelved workspace file for left hemisphere...")
root.withdraw()

try:
    
    exec(open('./RestoreShelvedWorkspaceScript.py').read())

except:
    
    print('Error when restoring shelved workspace.  Will attempt to continue')
    
CellFluorTraces_Frame_LH = CellFluorTraces_Frame
SessionID_LH = File[0:19]

root = tk.Tk()

RestoreFilePath = askopenfilename(initialdir=DefaultWorkspaceDir, 
                                  title="Select a shelved workspace file for right hemisphere...")
root.withdraw()

try:
    
    exec(open('./RestoreShelvedWorkspaceScript.py').read())

except:
    
    print('Error when restoring shelved workspace.  Will attempt to continue')
    
CellFluorTraces_Frame_RH = CellFluorTraces_Frame

(NumColumns,) = CellFluorTraces_Frame_RH.columns.values.shape

CellFluorTraces_Frame_LH = \
    CTDFL.TraceSampler(CellFluorTraces_Frame_LH, NumColumns - 1)

CellFluorTraces_Frame_Joined = \
    CTDFL.JoinDataFramesAcrossColumns(CellFluorTraces_Frame_LH, 
                                      CellFluorTraces_Frame_RH)
    
# Determine parent directory and filename from complete path.
Drive, Path_and_file = os.path.splitdrive(RestoreFilePath)
Path, File = os.path.split(Path_and_file)

# Extract session id from filename
SessionID = File[0:19]

# Construct a  workspace filename to use as save default
SaveRootDir = DefaultWorkspaceDir
SaveFilename =  SessionID + '_new_unique_400ms_SamFiltered_LH+RH'
DefaultSavePath = SaveRootDir + os.path.sep + SessionID 

CandidateSavePath = DefaultSavePath + os.path.sep + SaveFilename
result = askyesno(title="Save to default?", 
                  message="Save to default location:\n" + CandidateSavePath)

if result == True:
    
    # Check if default directory exists, if not, make it.
    if not os.path.exists(DefaultSavePath):
    
        os.makedirs(DefaultSavePath)
    
    # Set save path
    SavePath = CandidateSavePath

else:
    
    # Prompt user to enter the desired save path
    SavePath = asksaveasfilename(title='Set workspace save path',
                                 initialdir=DefaultSavePath,
                                 initialfile=SaveFilename)

# End file dialogs
root.withdraw()

# Run sliding window analysis using above-defined subroutine.
SWAF.SlidingWindowAnalysisFunc(BehavDict, CellFluorTraces_Frame_Joined, 
                          SavePath, ParamsDict)
# remember to  change session ID in files
SaveFilename =  SessionID_LH + '_new_unique_400ms_SamFiltered_SubsampledLH'
DefaultSavePath = SaveRootDir + os.path.sep + SessionID_LH 
SavePath2 = DefaultSavePath + os.path.sep + SaveFilename

SWAF.SlidingWindowAnalysisFunc(BehavDict, CellFluorTraces_Frame_LH, 
                          SavePath2, ParamsDict)
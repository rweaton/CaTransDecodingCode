#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 15:27:19 2019

@author: thugwithyoyo
"""

### Appears to work properly now. Shelve would not save the variable
### named: 'path'.  Fixed: 2019-12-02
#########################################################################
### !!!!! This program corrupts shelve files; fix before using !!!!!!####
#########################################################################


from PeriEventTraceFuncLib import *
import numpy as np
from CaTraceNormalizer import *
from collections import defaultdict
import os

import CalciumImagingFluorProcessing as CIFP

import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()


def RepairCellFluorTracesLabelsAddCentroidFrameFunc(RestoreFilePath):

    # Determine parent directory and filename from complete path.
    Drive, Path_and_file = os.path.splitdrive(RestoreFilePath)
    Path, File = os.path.split(Path_and_file)
    
    # Note: all variables below should be part of the restored workspace
    
    # Regenerate specified workspace
    exec(open('./RestoreShelvedWorkspaceScript.py').read())
    
    # Reproduce the calcium fluorescence frame; this time with correct labels!
    CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)
    
    # Generate dataframe of centroids 
    CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)
    
    # Save over previous workspace with the new and repaired dataframes.
    SavePath = RestoreFilePath[0:-4]
    
    #####################################
    ########   Start shelving   #########
    #####################################
    
    # Assemble save path
    # get root directory of save path from path to calcium data
    # SavePath is an argument above.  The following script requires it 
    exec(open('./ShelveWorkspaceScript.py').read())
    
RepairCellFluorTracesLabelsAddCentroidFrameFunc(RestoreFilePath)
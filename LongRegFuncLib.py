#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  5 19:11:22 2020

@author: thugwithyoyo
"""
#from PeriEventTraceFuncLib import *

import numpy as np
import pandas as pd
import PeriEventTraceFuncLib as PETFL

import os
import shelve
from collections import defaultdict

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename
from tkinter.messagebox import askyesno

#### Longitudinal Registration function library #### 

def BuildLongRegDataframeFromCSV(**kwargs):
    
    # Acquire path of workspace to load.
    # Set defualt parent directories
    DataRootDir = '/home/thugwithyoyo/CaTransDecoding/LongitudinalRegistrationLists'

    # Start file dialogs 
    root = tk.Tk()

    if kwargs['path_to_long_cellIDs'] is not None:
        
        lr_path = kwargs['path_to_long_cellIDs']
        
    else:
        
        # Prompt user to navigate to, and select, the longitudinal cell identities file
        lr_path = askopenfilename(title='Select longitudinal cell identities file',
                                      filetypes=[("csv files", "*.csv")],
                                      initialdir=DataRootDir)

    # End file dialogs
    root.withdraw()
    
    # Asssemble longitudinal registration information into pandas datatrame
    lr_frame = pd.read_csv(lr_path)
    
    return lr_frame

def BuildSessionDatasets(**kwargs):
   
    # Acquire path of workspace to load.
    # Set defualt parent directories
    DataRootDir = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData'
    
    # Start file dialogs 
    root = tk.Tk()
    
    if kwargs['path_to_behav_file'] is not None:
        
        # If a path is included as an argument, set it to the path variable
        PathToBehavFile = kwargs['path_to_behav_file']

    else:
        
        # Prompt user to navigate to, and select, the session behavior file
        PathToBehavFile = askopenfilename(title='Select session behavior file',
                                          filetypes=[("json files","*.json"), 
                                                     ("csv files", "*.csv")],
                                          initialdir=DataRootDir)
        
    if kwargs['path_to_fluor_file'] is not None:

        PathToFluorFile = kwargs['path_to_fluor_file']
        
    else:
        
        # Prompt user to navigate to, and select, the session imaging file
        PathToFluorFile = askopenfilename(title='Select session fluorescence record file',
                                          filetypes=[("json files","*.json"),
                                                     ("csv files", "*.csv")],
                                          initialdir=DataRootDir)

    # Assemble the dataframe of cell fluorescence traces.
    CellFluorTraces_Frame = PETFL.CellFluorTraces_FrameGen(PathToFluorFile)
    
    # Assemble the dictionary of event occurence timestamps
    BehavDict = PETFL.BehavDictGen(PathToBehavFile)

    # End file dialogs
    root.withdraw()

    # Determine parent directory and filename from complete path.
    Drive, Path_and_file = os.path.splitdrive(PathToBehavFile)
    Path, File = os.path.split(Path_and_file)
    
    # Extract session id from filename
    SessionID = File[0:19]
    
    return CellFluorTraces_Frame, BehavDict, SessionID


def ElementLocator(TargetsArray, TestArray):
    
    # This program finds the location of TargetsArray elements within the 
    # array TestArray, provided that they exist.  More specifically, for each
    # element in TargetsArray, the routine searches for a match in TestArray
    # and returns the location of that match (within TestArray) by its array 
    # index. The program assumes that elements of TestArray are unique; if 
    # multiple matches are found (1-to-n mapping) only the index of the first 
    # match is retained.
    
    # Count number of elements in TargetsArray
    (NumTargetElements,) = TargetsArray.shape
    
    # Count number of elements in  TestArray
    (NumTestElements,) = TestArray.shape
    
    # Generate sequential list of indices to count locations of elements in
    # TestArray
    TestArrayIndices = np.arange(0, NumTestElements)
    
    # Initialize the array to contain index locations to return
    TargetsMap = np.empty((NumTargetElements,), dtype=int)
    
    # Iterate over elements within TargetsArray
    for t in np.arange(0, NumTargetElements):
        
        # For each element, find matches within TestArray
        Filt = (TestArray == TargetsArray[t])
        
        # Write the index location of first match to the map output variable
        TargetsMap[t] = TestArrayIndices[Filt][0]
        
    return TargetsMap

def CommonCellsExtractor(lr_frame, Session1, Session2):
    
    # Generate filters to extract session subsets
    S1_filt = (lr_frame['local_cellset_file'] == Session1)
    S2_filt = (lr_frame['local_cellset_file'] == Session2)
    
    # Obtain cell lists by global index from the two sessions
    S1_global_cellset = lr_frame['global_cell_index'][S1_filt]
    S2_global_cellset = lr_frame['global_cell_index'][S2_filt]
    
    # Obtain cell lists by local cell name from  the two sessions
    S1_local_cellset = lr_frame['local_cell_name'][S1_filt]
    S2_local_cellset = lr_frame['local_cell_name'][S2_filt]
    
    # Determine particular cells that are present in both lists
    common_global_cellset = np.intersect1d(S1_global_cellset.values, 
                                           S2_global_cellset.values)
    
    # Sort the unique values of the result in ascending order
    common_global_cellset = np.sort(common_global_cellset)
    
    # Generate a map to cells, that are common to both sessions, within the
    # session 1 lists
    S1_cells_map = ElementLocator(common_global_cellset, 
                                    S1_global_cellset.values)
    
    # Generate a map to cells, that are common to both sessions, within the
    # session 2 lists
    S2_cells_map = ElementLocator(common_global_cellset,
                                    S2_global_cellset.values)
    
    return {
            'S1_cells_map' : S1_cells_map,
            'S2_cells_map' : S2_cells_map,
            'S1_global_cellset' : S1_global_cellset.values[S1_cells_map],
            'S2_global_cellset' : S2_global_cellset.values[S2_cells_map],
            'S1_local_cellset' : S1_local_cellset.values[S1_cells_map], 
            'S2_local_cellset' : S2_local_cellset.values[S2_cells_map]
            }
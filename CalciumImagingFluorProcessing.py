#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 07:11:38 2019

@author: thugwithyoyo
"""

import numpy as np
import os
import json
import pandas as pd

#PathToFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_001_Left_CellTraces.csv'

def IDPS_CellFluorTracesParser(PathToFile):

  # Aquire path of the directory that contains this script.  This is so 
  # the routine can navigate back after changing directories if need be.  
  ScriptDir = os.getcwd()
  
  
  CellFluorTraces_Frame = pd.read_csv(PathToFile, header=0, skiprows=[1])
  
  # Generate an index list for later referencing
  (nSamples, nColumns) = CellFluorTraces_Frame.shape
  
  # Generate the an array of indices for referencing
  Indices = np.arange(0, nSamples)
  
  # Properly name the first column; the list of time stamps.
  CellFluorTraces_Frame.columns.values[0] = 'Timestamps'
  
  # The list of column names.  The first column contains timestamps.  All
  # columns that follow are fluorescence traces of cells.
  ColumnNamesList = list(CellFluorTraces_Frame.columns.values)
  
  return CellFluorTraces_Frame

def FrameFromJSON(PathToFluorFile):
    
    # This routine opens a json file (converted from a .mat file) that contains
    # the calcium trace data and parses it into a pandas dataframe for sub-
    # sequent processing.
    
    # Load data from json file.
    with open(PathToFluorFile) as f:
        FluorDataList = json.load(f)
    
    # Count number of cells in file.    
    nCells = len(FluorDataList)
    #Overrode above defn for  2018-12-14 file.
    #nCells = 95
    
    # Determine number of samples in each trace.
    nSamplesPerTrace = len(FluorDataList[0]['trace_ts'])
        
    # Assemble list of column names and data array of column vectors
    # Initialize data array.
    Data = np.empty((nSamplesPerTrace, nCells + 1), dtype=float)
    Data[:,0] = FluorDataList[0]['trace_ts']
    
    # Initialize list to contain column names.
    ColumnNames = [""]*(nCells + 1)
    ColumnNames[0] = "Timestamps"
    CellNameWidth = len(str(nCells))

    # Populate the data array column-by-column and build corresponding name
    # of column. Except for the first column, each column after contains the
    # calcium trace of a cell.
    for i in np.arange(0,nCells):
        
        # Build cell id name.
        #ColumnNames[i + 1] = 'C' + str(i).zfill(CellNameWidth)
        
        # Acquire cell id names
        ColumnNames[i + 1] = FluorDataList[i]['name']
        
        # Populate column with calcium fluorescence trace.
        Data[:, i + 1] = FluorDataList[i]['trace'][0:nSamplesPerTrace]
       
        # This was output for debugging.
        #print(i)
        
    # Write data and column names array in to a pandas dataframe.
    return pd.DataFrame(data=Data, columns=ColumnNames)
    
def CellCentroidsFromJSON(PathToFluorFile):
    
    # This routine opens a json file (converted from a .mat file) that contains
    # the centroid coordiante data of each cell and parses it into a pandas 
    # dataframe for subsequent processing.
    
    # Load data from json file.
    with open(PathToFluorFile) as f:
        FluorDataList = json.load(f)
    
    # Count number of cells in file.    
    nCells = len(FluorDataList)
    
    # Assemble list of column names and data array of column vectors
    # Initialize data array.
    Data = np.empty((2, nCells), dtype=float)
    
    # Initialize list to contain index names
    IndexNames = np.array(['x_coord', 'y_coord'])

    # Initialize list to contain column names.
    ColumnNames = [""]*(nCells)

    # Populate the data array column-by-column and build corresponding name
    # of column. Except for the first column, each column after contains the
    # centroid location of a cell in the FOV.
    for i in np.arange(0,nCells):
        
        # Build cell id name.
        #ColumnNames[i] = 'C' + str(i).zfill(CellNameWidth)
        
        # Acquire cell id name
        ColumnNames[i] = FluorDataList[i]['name']
        
        # Populate column with calcium fluorescence trace.
        Data[:, i] = np.array(FluorDataList[i]['centroid'])
        
    return pd.DataFrame(data=Data, columns=ColumnNames, index=IndexNames)
   

def PeriEventTraceExtractor(CellFluorTraces_Frame, RefEvents, BoundaryWindow):

  # The list of column names.  The first column contains timestamps.  All
  # columns that follow are fluorescence traces of cells.
  ColumnNamesList = list(CellFluorTraces_Frame.columns.values)  
    
  # Generate an index list for later referencing
  (nSamples, nColumns) = CellFluorTraces_Frame.shape
  
  # Cell columns indices.  The first Column contains timestamps
  CellColumnsIndices = np.arange(1, nColumns)
  
  # Extract array of timestamps.
  Timestamps = CellFluorTraces_Frame['Timestamps'].values
  
  CalciumTimeInc = np.round(np.mean(np.diff(Timestamps)), decimals=8)
  
  
#  # Check that timestamps are evenly sampled
#  Incs = np.diff(Timestamps)
#  
#  if (np.sum(Incs != Incs[0]) != 0):
#      
#      print('Traces were not sampled evenly.  Exiting now.')
#      break
#  
#  else:
#      
#      continue
  #### Incorrect approach #############    
  # Generate row extraction filter.  Note: Though recordings from which  
  # event timestamps and calcium imaging began at the same time (synced),
  # the two sources were recorded using DIFFERENT CLOCKS at DIFFERENT 
  # RATES.  Accordingly, event timestamps must be "rounded" to match the 
  # sampling rate of calcium imaging (20 Hz?).
  #RefEvents = np.round(RefEvents/CalciumTimeInc)*CalciumTimeInc
  
  # Eliminate possible "jitter" in sampling of calcium activity signal
  #Timestamps = np.round(Timestamps/CalciumTimeInc)*CalciumTimeInc

  # Determine the number of samples that are contained in the extraction window.
  nSamplesPerTrace = int(np.ceil((BoundaryWindow[1] - BoundaryWindow[0])/
                              CalciumTimeInc))  
  
  # Initialize array to contain (flattened) array of traces
  nCells = nColumns - 1

  PeriEventActivity_flat = np.empty((RefEvents.size, nCells*nSamplesPerTrace))
  #PeriEventActivity = np.empty((RefEvents.size, nSamplesPerTrace, nCells))
  
  for RefEventIndex in np.arange(0, RefEvents.size):

      # Compute current window bounds
      CurrentWindow = BoundaryWindow + RefEvents[RefEventIndex]      
      
      # Generate a boolean filter for extracting peri-event calcium activity. 
      # NOTE: left window bound in inclusive, right window bound is exclusive.
      RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
          ( Timestamps < CurrentWindow[1] )
      
#      if np.sum(RowFilt) == (nSamplesPerTrace - 1):
#          
#          print("Short snippet detected. Appending to end one more element...")
#          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
#                ( Timestamps < (CurrentWindow[1] + CalciumTimeInc))
#                
#      elif np.sum(RowFilt) == (nSamplesPerTrace + 1):         
#
#          print("Long snippet detected.  Removing last element...")
#          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
#                ( Timestamps < (CurrentWindow[1] - CalciumTimeInc))
                
      while np.sum(RowFilt) < nSamplesPerTrace:
          
          print("Short snippet detected. Appending to end one more element...")
          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
                ( Timestamps < (CurrentWindow[1] + 0.5*CalciumTimeInc))
                
      while np.sum(RowFilt) > nSamplesPerTrace:
          
          print("Long snippet detected.  Removing last element...")
          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
                ( Timestamps < (CurrentWindow[1] - 0.5*CalciumTimeInc))
                
      # Populate peri-event array from dataframe.
      for CellIndex  in np.arange(1, nCells + 1):

          ColumnsSlice = slice((CellIndex - 1)*nSamplesPerTrace, CellIndex*nSamplesPerTrace, 1)
          PeriEventActivity_flat[RefEventIndex, ColumnsSlice] = \
              CellFluorTraces_Frame[ColumnNamesList[CellIndex]][RowFilt]

  return PeriEventActivity_flat
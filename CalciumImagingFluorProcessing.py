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
    
    # Load data from json file.
    with open(PathToFluorFile) as f:
        FluorDataList = json.load(f)
        
    nCells = len(FluorDataList)
    #Overrode above defn for  2018-12-14 file.
    #nCells = 95
    
    nSamplesPerTrace = len(FluorDataList[0]['trace_ts'])
        
    # Assemble list of column names and data array of column vectors
    # Initialize arrays
    Data = np.empty((nSamplesPerTrace, nCells + 1), dtype=float)
    Data[:,0] = FluorDataList[0]['trace_ts']
    
    ColumnNames = [""]*(nCells + 1)
    ColumnNames[0] = "Timestamps"
    CellNameWidth = len(str(nCells))

    for i in np.arange(0,nCells):
        
        ColumnNames[i + 1] = 'C' + str(i).zfill(CellNameWidth)
        Data[:, i + 1] = FluorDataList[i]['trace'][0:nSamplesPerTrace]
        print(i)
    return pd.DataFrame(data=Data, columns=ColumnNames)
    

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
      
      if np.sum(RowFilt) == (nSamplesPerTrace - 1):
          
          print("Short snippet detected. Appending to end one more element...")
          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
                ( Timestamps < CurrentWindow[1] + CalciumTimeInc)
                
      elif np.sum(RowFilt) == (nSamplesPerTrace + 1):         

          print("Long snippet detected.  Removing last element...")
          RowFilt = ( Timestamps >= CurrentWindow[0] ) & \
                ( Timestamps < CurrentWindow[1] - CalciumTimeInc)
                
      # Populate peri-event array from dataframe.
      for CellIndex  in np.arange(1, nCells + 1):

          ColumnsSlice = slice((CellIndex - 1)*nSamplesPerTrace, CellIndex*nSamplesPerTrace, 1)
          PeriEventActivity_flat[RefEventIndex, ColumnsSlice] = \
              CellFluorTraces_Frame[ColumnNamesList[CellIndex]][RowFilt]

  return PeriEventActivity_flat
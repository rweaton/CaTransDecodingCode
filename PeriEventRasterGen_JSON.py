#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 18:32:50 2019

@author: thugwithyoyo
"""
import json
import numpy as np
from hwfunclib import *


#PathToDataFile = '/media/thugwithyoyo/5477956AC6C57C61/mat data/2018-12-05-12-10-05_Struct.json'
PathToDataFile = '/media/thugwithyoyo/DA3F-23EC/mat data/2018-12-05-12-10-05_Struct.json'

RequiredArraysList = np.array([
        'Timestamp', 'Frame_Number', 'X1_1_pix', 'Y1_1_pix',
        'X2_1_pix', 'Y2_1_pix', 'EV1_1', 'EV1_3', 'EV1_7', 'EV1_9'])


# Set bin resolution to 1 ms for spike raster plot
Resolution = 0.005  # units in seconds

BoundaryWindow = (-5, 2)  # units in seconds

with open(PathToDataFile) as f:
  DataDict = json.load(f)
  
# Extract dict values as numpy arrays for subsequent processing
ArrayLength = len(DataDict['ReachData'])
ArrayIndices = np.arange(0, ArrayLength)
ArrayNames = list(DataDict['ReachData'][0].keys())

for ArrayName in ArrayNames:
    exec("%s = %s" % (ArrayName, 'np.zeros([ArrayLength])'))
    
for ArrayName in ArrayNames:
    
    for i in ArrayIndices:
        
        exec("RHS = DataDict['ReachData'][i]['%s']" % ArrayName)
        exec("%s[i] = %s" % (ArrayName, RHS))
        
nBuiltArrays = np.sum(np.isin(np.array(dir()), RequiredArraysList))
    
if nBuiltArrays == RequiredArraysList.size:
    print('All required arrays from ReachData have been generated.')
else:
    print ('There are required arrays that were not generated.')
    exit()
    
# Extract reach events
EV1_1_pos_events = Timestamp[np.hstack([False, np.diff(EV1_1) == 1])]
EV1_1_neg_events = Timestamp[np.hstack([False, np.diff(EV1_1) == -1])]

EV1_9_pos_events = Timestamp[np.hstack([False, np.diff(EV1_9) == 1])]
EV1_9_neg_events = Timestamp[np.hstack([False, np.diff(EV1_9) == -1])]

# Plot rasters
nROIs = len(DataDict['CalciumData'])
#nROIs= 125

nRows = 5
nCols = 5

for i in np.arange(0, nROIs):
#for i in np.arange(0, 10):
    
    if i % (nRows*nCols) == 0:
        
        # Set up subplot grid of size [nRows X nCols]
        fig, axs = plt.subplots(nrows=nRows, ncols=nCols)
        RowIndex = 0
        ColumnIndex = 0

    # Write group title to figure
    fig.suptitle('Peri-event spike rasters (bin-width = ' 
                 + str(np.round(Resolution, decimals=4)) + ' sec.)', fontsize=16)
    
    # Extract peri-event spike times
    PEtimes = PeriEventTimeStampExtractor(EV1_9_pos_events, 
                                          DataDict['CalciumData'][i]['events'], 
                                          BoundaryWindow)
        
    # Determine location of current subplot in grid
    j = i - (nRows*nCols)*int(np.floor(i/(nRows*nCols)))
    RowIndex = int(np.floor(j/nRows))
    
    if (i % nRows == 0):
        ColumnIndex = 0
        
    #RowIndex = EventIndices[Events == e][0]
    #ColumnIndex = NeuronIndices[Neurons == n][0]
        
    # Compile current raster
    RastDict = RasterMatrixCompiler(PEtimes, BoundaryWindow, Resolution)

    # Plot current raster
    CurrentAxs = axs[RowIndex, ColumnIndex]
    RasterPlot(RastDict, CurrentAxs, (RowIndex == nRows-1), 
                                     (ColumnIndex == 0))
    
    ColumnIndex += 1

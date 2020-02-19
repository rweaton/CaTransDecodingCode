#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 09:18:05 2020

@author: thugwithyoyo
"""
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

import os
import shelve

import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'

# Have user select workspace to load.  Workspace must be a shelve object '.dat'
# file generated from DecodingStabilityBatchProcessorFunc.py
DefaultWorkspaceDir = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing'

root = tk.Tk()

RestoreFilePath = askopenfilename(initialdir=DefaultWorkspaceDir, 
                                  title="Select a shelved workspace file...")
root.withdraw()

# Specify variables to load from  shelved workspace
VarsToLoad = ['DecodingStabilityOutput', 
              'DecodingSessionsList', 
              'LongSessionsList', 
#              'NumLatents', 
#              'NumShuffles', 
              'ParamsDict', 
              'SavePath']

#

# Save path of directory that contains this script
ScriptDir = os.getcwd()

(PathToFile, Filename) = os.path.split(RestoreFilePath)

# Change to directory that contains file to be loaded.
os.chdir(PathToFile)
 
# Open a shelf dictionary object that contains stored variables
my_shelf = shelve.open(os.path.splitext(Filename)[0])


#NumOfCommonCellsMat[i, j] = RunTestOnTrainedDict['NumOfCommonCells']
#TrainingSetPerfMat[i, j] = RunTestOnTrainedDict['TrainingSetPerf']
#TestingSetPerfMat[i, j] = RunTestOnTrainedDict['TestingSetPerf']
#OutcomeShufTrainingSetPerfMat[i, j, :] = RunTestOnTrainedDict['OutcomeShufTrainingSetPerfVec']
#OutcomeShufTestingSetPerfMat[i, j, :] = RunTestOnTrainedDict['OutcomeShufTestingSetPerfVec']
#CellShufTrainingSetPerfMat[i, j, :] = RunTestOnTrainedDict['CellShufTrainingSetPerfVec']
#CellShufTestingSetPerfMat[i, j, :] = RunTestOnTrainedDict['CellShufTestingSetPerfVec']   


# Iterate through list contents and write each to the global variables
# list.
for i in np.arange(0, len(VarsToLoad)):
    exec(VarsToLoad[i] + '= my_shelf[VarsToLoad[{}]]'.format(i))

# Close the shelve object
my_shelf.close()

# Return to script directory
os.chdir(ScriptDir)

###############################################################################
########                     BEGIN PLOT PROCESSING                   ##########       
###############################################################################
#### Compute durations (in units of days) separating sessions ####
# Initialize the days separation matrix
DaysSeparationMat = np.empty_like(DecodingStabilityOutput, dtype=float)

# Count number of sessions in list
(NumSessions,) = DecodingSessionsList.shape

# Generate index list for iterating through sessions
SessionsIndices = np.arange(0, NumSessions)

# Retrieve number of shuffles for matrix initialization
NumShuffles = ParamsDict['NumRepetitions']

# Calculate constant for seconds to days conversion
SecondsInADay =  24 * 3600.

# Initialize matrices (to be flattened) to plot
NumOfCommonCellsMat = np.empty((NumSessions, NumSessions), dtype=int)

TrainingSetCVPerfMat = np.empty((NumSessions, NumSessions), dtype=float)

TrainingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)

TestingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)

OutcomeShufTrainingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                  dtype=float)

OutcomeShufTestingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                  dtype=float)    

CellShufTrainingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                  dtype=float)

CellShufTestingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                  dtype=float)

# Iterate through all possible combinations of sessions list calculating 
# days separation between each pair of sessions
for i in SessionsIndices:
    
    for j in SessionsIndices:
        
        # Import the session time labels into datetime objects
        Sess1_dt = datetime.strptime(DecodingSessionsList[i], "%Y-%m-%d-%H-%M-%S")
        Sess2_dt = datetime.strptime(DecodingSessionsList[j], "%Y-%m-%d-%H-%M-%S")
        
        # Compute the duration of time between the two sessions by subtracting 
        # the second from the first.
        SessDiff_dt = Sess2_dt - Sess1_dt
        
        # Convert the separation duration into units of days and save to
        # appropriate location of storage matrix
        DaysSeparationMat[i, j] = \
            SessDiff_dt.days + SessDiff_dt.seconds/SecondsInADay
            
        NumOfCommonCellsMat[i, j] = DecodingStabilityOutput[i, j]['NumOfCommonCells']
        
        TrainingSetCVPerfMat[i, j] = DecodingStabilityOutput[i, j]['TrainingSetCVPerf']
        TrainingSetPerfMat[i, j] = DecodingStabilityOutput[i, j]['TrainingSetPerf']
        TestingSetPerfMat[i, j] = DecodingStabilityOutput[i, j]['TestingSetPerf']
        
        OutcomeShufTrainingSetPerfMat[i, j, :] = DecodingStabilityOutput[i, j]['OutcomeShufTrainingSetPerfVec']
        OutcomeShufTestingSetPerfMat[i, j, :] = DecodingStabilityOutput[i, j]['OutcomeShufTestingSetPerfVec']
        CellShufTrainingSetPerfMat[i, j, :] = DecodingStabilityOutput[i, j]['CellShufTrainingSetPerfVec']
        CellShufTestingSetPerfMat[i, j, :] = DecodingStabilityOutput[i, j]['CellShufTestingSetPerfVec'] 

# Compute y value proportions (test set perf. / training set perf.) and collapse
# 2D array over rows.
#yVals =  np.divide(DecodingStabilityDict['TestingSetPerfMat'], 
#                   DecodingStabilityDict['TrainingSetPerfMat']).flatten(order='C')

# Acquire example session indices
ExTrainSet = '2018-12-18-11-20-21'
ExTestSet = '2018-12-17-11-38-42'

ExTrainSetIndex = SessionsIndices[DecodingSessionsList == ExTrainSet]
ExTestSetIndex = SessionsIndices[DecodingSessionsList == ExTestSet]

# Extract example session data point
xVal_ex = DaysSeparationMat[ExTrainSetIndex, ExTestSetIndex]
yVal_ex = np.divide(TestingSetPerfMat[ExTrainSetIndex, ExTestSetIndex], 
                    TrainingSetCVPerfMat[ExTrainSetIndex, ExTestSetIndex])

# Set example session datapoint to nan in compliment array so it can be 
# plotted separately
TestingSetPerfMat[ExTrainSetIndex, ExTestSetIndex] = np.nan

# TFB CV performance over CV performance for within session datapoints 
# to make co-authors happy... TFB=Total Fucking Bullshit!!!
for i in SessionsIndices:
    
    TestingSetPerfMat[i, i] = TrainingSetCVPerfMat[i, i]
            
yVals =  np.divide(TestingSetPerfMat, 
                   TrainingSetCVPerfMat).flatten(order='C')

#CellShufRatio = \
#    np.divide(DecodingStabilityDict['CellShufTestingSetPerfMat'],  
#              np.repeat(DecodingStabilityDict['TrainingSetPerfMat'][:, :, np.newaxis], 
#                        NumShuffles, axis=2))

CellShufRatio = \
    np.divide(CellShufTestingSetPerfMat,  
              np.repeat(TrainingSetCVPerfMat[:, :, np.newaxis], 
                        NumShuffles, axis=2))
    
yVals_ShufCells = \
    np.mean(CellShufRatio, axis=2).flatten(order='C')
    
yVals_ShufCells_SE = \
    (1./np.sqrt(NumShuffles))*np.std(CellShufRatio, axis=2, ddof=1).flatten(order='C')
    
#TrialsShufRatio = \
#    np.divide(DecodingStabilityDict['OutcomeShufTestingSetPerfMat'],  
#              np.repeat(DecodingStabilityDict['TrainingSetPerfMat'][:, :, np.newaxis], 
#                        NumShuffles, axis=2))
    
TrialsShufRatio = \
    np.divide(OutcomeShufTestingSetPerfMat,  
              np.repeat(TrainingSetCVPerfMat[:, :, np.newaxis], 
                        NumShuffles, axis=2))
    
yVals_ShufTrialsModel = \
    np.mean(TrialsShufRatio, axis=2).flatten(order='C')

yVals_ShufTrialsModel_SE = \
    (1./np.sqrt(NumShuffles))*np.std(TrialsShufRatio, axis=2, ddof=1).flatten(order='C')

# collapse over rows the days separation array.  The flattened array will define
# the x axis values of datapoints in the scatterplot.
xVals = DaysSeparationMat.flatten(order='C')

figure, axs = plt.subplots(nrows=1, ncols=1)

axs.scatter(xVals, yVals, marker='o', label='observed outcomes', color='black')

# Plot example session datapoint
axs.scatter(xVal_ex, yVal_ex, marker='o', label='observed outcome (example)', 
            color='red')

#axs.errorbar(xVals, yVals_ShufCells, yerr=yVals_ShufCells_SE, marker='o', fmt='.',
#             label='shuffled cells', color='gray')


axs.errorbar(xVals, yVals_ShufTrialsModel, yerr=yVals_ShufTrialsModel_SE,
             marker='o', fmt='.', 
             label='shuffled outcomes model', color='lightblue')

axs.set_xlabel('intersession time (days)')
axs.set_ylabel('performance ratio')

axs.legend()

StartDate = LongSessionsList[0][0:10]
EndDate = LongSessionsList[-1][0:10]
axs.set_title(StartDate + ' through ' + EndDate + ', (' + str(NumSessions) + 
          r' sessions), N$_{latents}$ = ' + str(ParamsDict['NumLatents']))

figure.suptitle('Trace decoding stability')
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 07:22:37 2020

@author: thugwithyoyo
"""

import pandas  as pd
import numpy as np

# Table generator
#            'NumOfCommonCells': NumOfCommonCells,
#            'TrainingSetPerf': TestOnTrainedPerfDict['TrainingSetPerf'],
#            'TestingSetPerf': TestOnTrainedPerfDict['TestingSetPerf'],
#            'OutcomeShufTrainingSetPerfVec': OutcomeShufTrainingSetPerfVec,
#            'OutcomeShufTestingSetPerfVec': OutcomeShufTestingSetPerfVec,
#            'CellShufTrainingSetPerfVec': CellShufTrainingSetPerfVec,
#            'CellShufTestingSetPerfVec': CellShufTestingSetPerfVec  
#            
#VarsToSave = ['DecodingStabilityOutput', 
#              'DecodingSessionsList', 
#              'LongSessionsList', 
##              'NumLatents', 
##              'NumShuffles', 
#              'ParamsDict', 
#              'SavePath']

# Begin script
def GenerateBatchDecodingStabilityTable(DecodingStabilityOutput, 
                                        DecodingSessionsList, ParamsDict):
    
    (NumSessions,) = DecodingSessionsList.shape
    
    SessionsIndices = np.arange(0, NumSessions)
    
    #DecodingStabilityOutput = DecodingStabilityDict
    
    # Initialize dataframe to contain table
    DecodingStabilityTable = pd.DataFrame(columns=[
            'NumLatents', 'TrainingSet', 'TestingSet', 'NumOfCommonCells',
            'TrainingSetCVPerf', 'TrainingSetPerf', 'TestingSetPerf', 
            'TrainedOnTestSetCVPerf', 'TrainedOnTestSetPerf',
            'OutcomeShufTrainingSetPerf_Mean', 'OutcomeShufTrainingSetPerf_SEM',
            'OutcomeShufTestingSetPerf_Mean', 'OutcomeShufTestingSetPerf_SEM',
            'CellShufTrainingSetPerf_Mean', 'CellShufTrainingSetPerf_SEM',
            'CellShufTestingSetPerf_Mean', 'CellShufTestingSetPerf_SEM'
                                                   ],
                                          index=np.arange(0, NumSessions**2))
    
    # Begin script
    # Extract necessary parameters
    NumShuffles = ParamsDict['NumRepetitions']
    NumLatents = ParamsDict['NumLatents']
    
    IndexCounter = 0
    
    for i in SessionsIndices:
        
        for j in SessionsIndices:
            
            # Write parameters into table row for the current interation
            DecodingStabilityTable['NumLatents'][IndexCounter] = NumLatents
            
            DecodingStabilityTable['TrainingSet'][IndexCounter] = \
                DecodingSessionsList[i]
            
            DecodingStabilityTable['TestingSet'][IndexCounter] = \
                DecodingSessionsList[j]
            
            # Include number of common cells in the training/testing pair
            DecodingStabilityTable['NumOfCommonCells'][IndexCounter] = \
                DecodingStabilityOutput[i,j]['NumOfCommonCells']
            
            # Write performance into table row for the current iteration
            DecodingStabilityTable['TrainingSetCVPerf'][IndexCounter] = \
                DecodingStabilityOutput[i,j]['TrainingSetCVPerf']
            
            DecodingStabilityTable['TrainingSetPerf'][IndexCounter] = \
                DecodingStabilityOutput[i,j]['TrainingSetPerf']
                
            DecodingStabilityTable['TestingSetPerf'][IndexCounter] = \
                DecodingStabilityOutput[i,j]['TestingSetPerf']
    
            # NOTE the transposition of indices here!!!
            DecodingStabilityTable['TrainedOnTestSetCVPerf'][IndexCounter] = \
                DecodingStabilityOutput[j,i]['TrainingSetCVPerf']
                
            DecodingStabilityTable['TrainedOnTestSetPerf'][IndexCounter] = \
                DecodingStabilityOutput[j,i]['TrainingSetPerf']
    
            # Write performance of shuffling runs for the current iteration
            # Outcome shuffled, training set
            CurrentMean = np.mean(
                    DecodingStabilityOutput[i,j]['OutcomeShufTrainingSetPerfVec'])
            
            CurrentSE = np.std(
                    DecodingStabilityOutput[i,j]['OutcomeShufTrainingSetPerfVec'])/np.sqrt(NumShuffles)
            
            DecodingStabilityTable['OutcomeShufTrainingSetPerf_Mean'][IndexCounter] = \
               CurrentMean
            
            DecodingStabilityTable['OutcomeShufTrainingSetPerf_SEM'][IndexCounter] = \
               CurrentSE
               
            # Outcome shuffled, testing set
            CurrentMean = np.mean(
                    DecodingStabilityOutput[i,j]['OutcomeShufTestingSetPerfVec'])
            
            CurrentSE = np.std(
                    DecodingStabilityOutput[i,j]['OutcomeShufTestingSetPerfVec'])/np.sqrt(NumShuffles)
            
            DecodingStabilityTable['OutcomeShufTestingSetPerf_Mean'][IndexCounter] = \
               CurrentMean
            
            DecodingStabilityTable['OutcomeShufTestingSetPerf_SEM'][IndexCounter] = \
               CurrentSE
            
            # Cell shuffled, training set
            CurrentMean = np.mean(
                    DecodingStabilityOutput[i,j]['CellShufTrainingSetPerfVec'])
            
            CurrentSE = np.std(
                    DecodingStabilityOutput[i,j]['CellShufTrainingSetPerfVec'])/np.sqrt(NumShuffles)
            
            DecodingStabilityTable['CellShufTrainingSetPerf_Mean'][IndexCounter] = \
               CurrentMean
            
            DecodingStabilityTable['CellShufTrainingSetPerf_SEM'][IndexCounter] = \
               CurrentSE        
            
            # Cell shuffled, testing set
            CurrentMean = np.mean(
                    DecodingStabilityOutput[i,j]['CellShufTestingSetPerfVec'])
            
            CurrentSE = np.std(
                    DecodingStabilityOutput[i,j]['CellShufTestingSetPerfVec'])/np.sqrt(NumShuffles)
            
            DecodingStabilityTable['CellShufTestingSetPerf_Mean'][IndexCounter] = \
               CurrentMean
            
            DecodingStabilityTable['CellShufTestingSetPerf_SEM'][IndexCounter] = \
               CurrentSE
               
            IndexCounter += 1
            
    return DecodingStabilityTable

DecodingStabilityTable = GenerateBatchDecodingStabilityTable(DecodingStabilityOutput, 
                                        DecodingSessionsList, ParamsDict)

SpreadsheetSaveName = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/DecSet_LeftHemRightArm_5latents_3.xlsx'
DecodingStabilityTable.to_excel(SpreadsheetSaveName)
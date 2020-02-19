#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 08:38:23 2020

@author: thugwithyoyo
"""

import numpy as np
import pandas as pd
import PeriEventTraceFuncLib as PETFL
import LongRegFuncLib as LRFL 

from collections import defaultdict

import os
from os import listdir
from os.path import isfile, join

import shelve

#from tkinter import *
#from tkinter import filedialog


def DecodingStabilityBatchProcessor(ParamsDict, *,
                                    already_parsed=False, **kwargs):
    
    global CellFluorTraces_Frame
#    CellFluorTraces_Frame = pd.DataFrame()
    
    global BehavDict
#    BehavDict = {}
    
    VarsToLoad = ['CellFluorTraces_Frame', 'BehavDict']
    
    # Aquire path of current directory
    ScriptDir = os.getcwd()
    
    # Build from file, a dataframe longitudinal registration information that
    # includes global cell idenities mapped to cell identities within local
    # sessions.
    lr_frame = LRFL.BuildLongRegDataframeFromCSV()
    
    # For this routine to function, each element in the longitudinal registration 
    # sessions list must have a corresponding element in decoding sessions 
    # list that was recorded on the same day from the same hemisphere.
    #if kwargs['long_reg_sessions_list'] is not None:
    if np.isin('long_reg_sessions_list', np.array(list(kwargs.keys()))):
        
        # Set the list of sessions to use for longitudinal registration mapping
        # to the value of the specified argument
        LongRegSessionsList = np.array(kwargs['long_reg_sessions_list'])
        
    else:
        
        # Generate list of N unique sessions included in the longitudinal cell
        # identities file (Sabrina's list).
        LongRegSessionsList = np.unique(lr_frame['local_cellset_file'].values)
    
    
    #if  kwargs['decoding_sessions_list'] is not None:
    if np.isin('decoding_sessions_list', np.array(list(kwargs.keys()))):
        
        # Set the list of sessions to use for decoding to the value of the 
        # specified argument.
        SessionsList = np.array(kwargs['decoding_sessions_list'])
        
    else:
        
        # Assign the longitudinal registration sessions list to the decoding 
        # sessions list
        SessionsList = LongRegSessionsList
    
    # Count number of sessions to process.  Note: the longitudinal registration
    # sessions list must be of equivalent length
    (NumSessions,) = SessionsList.shape
    
    # Specify path to directory that contains corresponding datafiles for each
    # entry in SessionsList
    #PathToDataDir = '/media/thugwithyoyo/SP UFD U3/LongRegBatchProcessing/Datafiles'
    PathToDataDir = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/DataFiles'
    
    # Specify path to directory that will contain shelved workspaces for each
    # session after parsing.
    #PathToParsedDataDir = '/media/thugwithyoyo/SP UFD U3/LongRegBatchProcessing/ParsedWorkspaces'
    PathToParsedDataDir = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/ParsedWorkspaces'

    # Pack into a dict the two lists that specify the names and assigned 
    # numerical values of the behavioral event types listed in ParamsDict.
    # This dictionary will be used in later calls to the PLS decoder.
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
    
    #### Involved file handling here NOT FINISHED YET ###
    # Examine directory containing datafiles for each of the sessions in list.
    # Compare with unique sessions list from above.

    
    #os.chdir(PathToDataDir)
    #onlyfiles = [f for f in listdir(PathToDataDir) if isfile(join(PathToDataDir, f))]
    
    # For each entry in the unique sessions list that is NOT in the directory
    # list, prompt user to navigate to the corresponding datafile.
    
    ####################################
    
    ########### BEGIN PARSING SESSION DATA INTO SHELVED WORKSPACES ############
    # Assume all corresponding datafiles are in data directory.
    # Iterate through each entry in  SessionsList, build path to the datafiles
    # for constructing behavioral and fluorescence data structures.  Shelve
    # datastructures to workspace files to be opened and used as necessary.
    if already_parsed is False:
    #if (~np.isin('already_parsed', np.array(list(kwargs.keys())))) or (kwargs['already_parsed'] == False):
        
        for Session in SessionsList:
            
            # Construct path to behavioral datafile
            BehavFile = Session + '_new_unique_B.json'
            PathToBehavFile = join(PathToDataDir, BehavFile)
            
            # Construct path to fluorescence datafile
            FluorFile = Session + '_new_unique_C.json'
            PathToFluorFile = join(PathToDataDir, FluorFile)
            
            # Build peri-event datasets for the current session
            CellFluorTraces_Frame, BehavDict, SessionName = \
                LRFL.BuildSessionDatasets(path_to_behav_file=PathToBehavFile,
                                          path_to_fluor_file=PathToFluorFile)
                
            # Shelve datasets generated above
            # Construct save path
            VarsToSave = ['CellFluorTraces_Frame', 'BehavDict', 'SessionName']
            SaveFile = Session + '_SessionData'
            SavePath = join(PathToParsedDataDir, SaveFile)
            
            # Initialize and open self object
            my_shelf = shelve.open(SavePath)
            
            # Iterate through the list of keys in the save list, writing each to shelf
            # object
            for key in VarsToSave:
                
                # Surround write command in try/except clause to bypass possible TypeError
                # issues.
                try:
                    
                    my_shelf[key] = globals()[key]
                    #my_shelf[key] = locals()[key]
                    
                except TypeError:
                    #
                    # __builtins__, my_shelf, and imported modules can not be shelved.
                    #
                    print('ERROR shelving: {0}'.format(key))
             
            # Close shelf object after variables have been written to file
            my_shelf.close()
            ###### End workspace save ###########
        
        ############### END PARSING LOOP ######################################
        
    ################ INITIALIZE STORAGE ARRAYS BLOCK ######################
    NumShuffles = ParamsDict['NumRepetitions']
#    NumLatents = ParamsDict['NumLatents']
    if kwargs['analysis_type'] == 'sliding_window':
        
#        SlidingWindowDictArray = np.empty((NumSessions, NumSessions), dtype=dict)
        
        Output = np.empty((NumSessions, NumSessions), dtype=np.ndarray)        
        
    else:
        
#        NumOfCommonCellsMat = np.empty((NumSessions, NumSessions), dtype=int)
#        
#        TrainingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)
#        
#        TestingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)
#        
#        OutcomeShufTrainingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
#                                          dtype=float)
#        
#        OutcomeShufTestingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
#                                          dtype=float)    
#        
#        CellShufTrainingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
#                                          dtype=float)
#        
#        CellShufTestingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
#                                          dtype=float)
        
        Output = np.empty((NumSessions, NumSessions), dtype=dict)
        
    #CellFluorTraces_Frame = pd.DataFrame()
    #BehavDict = {}
    
    ############ END STORAGE ARRAY INITIALIZATION BLOCK ###################
    
    
    ####################  BEGIN BATCH PROCESSING LOOP #####################
    SessionsListIndices = np.arange(0, NumSessions)
#    ShufflesListIndices = np.arange(0, NumShuffles)
    
    # Iterate over all possible combinations of training and testing sets
    for i in SessionsListIndices:

        for j in SessionsListIndices:
            
            # Have a heart; don't keep the poor user in suspense.
            print('Processing training set: ' + SessionsList[i] + 
                  ', testing set: ' + SessionsList[j])

            # Obtain list of common cells between current training set 
            # (index i) and testing set (index j).
            CommonCellsDict = LRFL.CommonCellsExtractor(lr_frame, 
                                        LongRegSessionsList[i], 
                                        LongRegSessionsList[j])
            
            # Retrieve parsed workspaces of current training and testing
            # datasets
    
            # Open TRAINING set parsed workspace
            RestoreFile = SessionsList[i] + '_SessionData.dat'
            RestoreFilePath = join(PathToParsedDataDir, RestoreFile)
            #exec(open('./RestoreShelvedWorkspaceScript_Specific.py').read())
                        
            (PathToFile, Filename) = os.path.split(RestoreFilePath)
            
            # Change to directory that contains file to be loaded.
            os.chdir(PathToFile)
             
            # Open a shelf dictionary object that contains stored variables
            my_shelf = shelve.open(os.path.splitext(Filename)[0])
            
            # Iterate through dictionary contents and write each to the global variables
            # list.
            for key in VarsToLoad:
                
                globals()[key]=my_shelf[key]
                #locals()[key]=my_shelf[key]
            
            # Close the shelve object
            my_shelf.close()
            
            # Return to script directory
            #os.chdir(ScriptDir)            
            
            #print('ParamsDict after loading workspace: ' + str(ParamsDict))
            
            # Assign the workspace fluorescence dataframe for TRAINING
            CellFluorTraces_Frame_train = CellFluorTraces_Frame
            
            # Assign behavioral data dict for TRAINING set
            BehavDict_train = BehavDict
            
            # Remove events that occured either too early, or too late, with 
            # respect to the timeframe of the calcium fluorescence record. 
            # That is, events that would otherwise cause extraction of 
            # surrounding trace snippets that would be incomplete (not of full 
            # length).
#            BehavDict_train = PETFL.EarlyAndLateBehavEventRemoval(BehavDict, 
#                                                      CellFluorTraces_Frame_train, 
#                                                      RefEventsDict, 
#                                                      ParamsDict['BoundaryWindow'])

            # Include only cells listed in common cells dictionary that are also
            # listed in the cell traces dataframe for training.
#            Mask1 = np.isin(CommonCellsDict['S1_local_cellset'], 
#                            CellFluorTraces_Frame_train.columns.values)
                        
            # Open TESTING set parsed workspace
            RestoreFile = SessionsList[j] + '_SessionData.dat'
            RestoreFilePath = join(PathToParsedDataDir, RestoreFile)
            #exec(open('./RestoreShelvedWorkspaceScript_Specific.py').read())
                        
            (PathToFile, Filename) = os.path.split(RestoreFilePath)
            
            # Change to directory that contains file to be loaded.
            os.chdir(PathToFile)
             
            # Open a shelf dictionary object that contains stored variables
            my_shelf = shelve.open(os.path.splitext(Filename)[0])
            
            # Iterate through dictionary contents and write each to the global variables
            # list.
            for key in VarsToLoad:
                
                globals()[key]=my_shelf[key]
                #locals()[key]=my_shelf[key]
            
            # Close the shelve object
            my_shelf.close()
            
            # Return to script directory
            os.chdir(ScriptDir)
            
            # Assign the workspace fluorescence dataframe for TESTING
            CellFluorTraces_Frame_test = CellFluorTraces_Frame
            
            # Assign behavioral data dict for TESTING set
            BehavDict_test = BehavDict
            
            # Remove events that occured either too early, or too late, with 
            # respect to the timeframe of the calcium fluorescence record.             
#            BehavDict_test = PETFL.EarlyAndLateBehavEventRemoval(BehavDict, 
#                                                      CellFluorTraces_Frame_test, 
#                                                      RefEventsDict, 
#                                                      ParamsDict['BoundaryWindow'])            
            
#            # Include only cells listed in common cells dictionary that are also
#            # listed in the cell traces dataframe for testing.
#            Mask2 = np.isin(CommonCellsDict['S2_local_cellset'], 
#                            CellFluorTraces_Frame_test.columns.values)
#
#            # Set intersect the two masks to obtain the final selection mask for
#            # the cell fluorescence sets for both the training and testing sessions
#            SelectFilt = (Mask1 & Mask2)
#            
#            # Count the number of common cells to be used for this combination
#            # of training and testing sessions.
#            NumOfCommonCellsMat[i, j] = np.sum(SelectFilt)
#            
#            # Select common cell traces from TRAINING set            
#            CellFluorTraces_Frame_train = CellFluorTraces_Frame_train[
#                    np.hstack([np.array(['Timestamps']), 
#                               CommonCellsDict['S1_local_cellset'][SelectFilt]])]
#            
#            # Select common cell traces from TESTING set
#            CellFluorTraces_Frame_test = CellFluorTraces_Frame[
#                    np.hstack([np.array(['Timestamps']), 
#                               CommonCellsDict['S2_local_cellset'][SelectFilt]])]
#            
#            # Run peri-event extraction and assemble the dictionary to be 
#            # used for PLS-regression training.
#            PeriEventExtractorDict_train = \
#                PETFL.PeriEventExtractor_Trace(BehavDict_train, 
#                                               CellFluorTraces_Frame_train, 
#                                               RefEventsDict, 
#                                               ParamsDict['BoundaryWindow'])
#                
#            # Run peri-event extraction and assemble the dictionary to be 
#            # used to test the PLS-regression model above.
#            PeriEventExtractorDict_test = \
#                PETFL.PeriEventExtractor_Trace(BehavDict_test, 
#                                               CellFluorTraces_Frame_test, 
#                                               RefEventsDict, 
#                                               ParamsDict['BoundaryWindow'])
#                
#            # Compute the PLS regression matrices from the TRAINING dataset.
#            Performance_train = \
#                PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
#                                             NumLatents)
#                
#            #TrainingSetPerfMat[i, j] = Performance_train['performance']
#            # Apply the train dataset to the current model trained above
#            TrainingSetPerfMat[i, j] = \
#                PETFL.ApplyPLSModel(Performance_train['B'], 
#                                    Performance_train['B_0'], 
#                                    PeriEventExtractorDict_train)['performance']            
#            
#            
#            # Apply the test dataset to the current model trained above
#            TestingSetPerfMat[i, j] = \
#                PETFL.ApplyPLSModel(Performance_train['B'], 
#                                    Performance_train['B_0'], 
#                                    PeriEventExtractorDict_test)['performance']
#                
#            for k in ShufflesListIndices:
#                
#                #### TRAINING SET ####
#                # Shuffle target outcomes of TRAINING set and retrain 
#                (nTotalTrials, nTotalFeatures) = \
#                    PeriEventExtractorDict_train['PEA_Array'].shape
#                   
#                # Generate a map to shuffle target outcomes
#                ShuffledIndices = np.arange(0, nTotalTrials)
#                np.random.shuffle(ShuffledIndices)
#                
#                # Shuffle target outcome correspondence to be passed to PLS trainer
#                PeriEventExtractorDict_train['TargetsVec'] = \
#                    PeriEventExtractorDict_train['TargetsVec'][ShuffledIndices]
#                    
#                # Train PLS model on shuffled outcomes data
#                ShufOutcomesPerformance_train = \
#                    PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
#                                                 NumLatents)
#                
#                # Record performance in storage matrix
#                OutcomeShufTrainingSetPerfMat[i, j, k] = \
#                    ShufOutcomesPerformance_train['performance']
#                
#                # Shuffle TRAINING set cell order and run model on the result
#                # Shuffle order of cell fluorescence traces in the dataframe
#                CellFluorTraces_Frame_train_shuf = \
#                    CellFluorTraces_Frame_train[
#                            np.hstack([np.array(['Timestamps']), 
#                                       np.random.choice(CommonCellsDict['S1_local_cellset'][SelectFilt], 
#                                                        size=NumOfCommonCellsMat[i, j], 
#                                                        replace=False)])]
#                
#                # Generate test data structure with shuffled cell order
#                # of the TRAINING set
#                PeriEventExtractorDict_train_shuf = \
#                    PETFL.PeriEventExtractor_Trace(BehavDict_train, 
#                                    CellFluorTraces_Frame_train_shuf, 
#                                    RefEventsDict, 
#                                    ParamsDict['BoundaryWindow'])
#    
#                # Run model on shuffled data structure and record performance
#                CellShufTrainingSetPerfMat[i, j, k] = \
#                    PETFL.ApplyPLSModel(Performance_train['B'], 
#                                        Performance_train['B_0'], 
#                                        PeriEventExtractorDict_train_shuf)['performance']
#                
#                #### TESTING SET ####
#                # Run TESTING set in shuffled-outcome model and record performance
#                OutcomeShufTestingSetPerfMat[i, j, k] = \
#                    PETFL.ApplyPLSModel(ShufOutcomesPerformance_train['B'], 
#                                        ShufOutcomesPerformance_train['B_0'], 
#                                        PeriEventExtractorDict_test)['performance']
#                                    
#                # Shuffle TESTING set cell order and run model on the result
#                # Shuffle order of cell fluorescence traces in the dataframe
#                CellFluorTraces_Frame_test_shuf = \
#                    CellFluorTraces_Frame_test[
#                            np.hstack([np.array(['Timestamps']), 
#                                       np.random.choice(CommonCellsDict['S2_local_cellset'][SelectFilt], 
#                                                        size=NumOfCommonCellsMat[i, j], 
#                                                        replace=False)])]
#                
#                # Generate test data structure with shuffled cell order
#                # of the TESTING set
#                PeriEventExtractorDict_test_shuf = \
#                    PETFL.PeriEventExtractor_Trace(BehavDict_test, 
#                                    CellFluorTraces_Frame_test_shuf, 
#                                    RefEventsDict, 
#                                    ParamsDict['BoundaryWindow'])
#    
#                # Run model on shuffled data structure and record performance
#                CellShufTestingSetPerfMat[i, j, k] = \
#                    PETFL.ApplyPLSModel(Performance_train['B'], 
#                                        Performance_train['B_0'], 
#                                        PeriEventExtractorDict_test_shuf)['performance']
        
            # Above code put into subroutines in LongRegFuncLib
            # Evaluate performance for the current traning set and test set 
            # combination
            if kwargs['analysis_type'] == 'sliding_window':
                
                #print('Batch Level Boundary window: ' + str(ParamsDict['BoundaryWindow']))
                Output[i, j] = \
                    LRFL.SlidingWindowRunTestOnTrained(ParamsDict, CommonCellsDict, 
                                                       CellFluorTraces_Frame_train, 
                                                       BehavDict_train, 
                                                       CellFluorTraces_Frame_test, 
                                                       BehavDict_test)
                
                print(Output[i, j])
                
            else:
                
                Output[i, j] = \
                    LRFL.RunTestOnTrained(ParamsDict, CommonCellsDict, 
                                          CellFluorTraces_Frame_train, BehavDict_train, 
                                          CellFluorTraces_Frame_test, BehavDict_test)
                        
    return Output
    
    ############################# END FUNCTION ################################
    
############################# START EXECUTION SCRIPT ##########################
    
LongSessionsList = np.array([
#                             '2018-12-04-10-31-21',
#                             '2018-12-10-11-37-56',
#                             '2018-12-11-10-53-04',
#                             '2018-12-14-11-01-41',
                             '2018-12-17-11-38-42',
                             '2018-12-18-11-20-21',
#                             '2018-12-28-10-50-15'
                             ])
    
DecodingSessionsList = np.array([
#                                 '2018-12-04-10-31-21',
#                                 '2018-12-10-11-37-56',
#                                 '2018-12-11-10-53-04',
#                                 '2018-12-14-11-01-41',
                                 '2018-12-17-11-38-42',
                                 '2018-12-18-11-20-21',
#                                 '2018-12-28-11-12-14'
                                 ])
    
# Define a dictionary to contain keyed access to sliding window analysis parameters
ParamsDict = defaultdict(dict)

# Specify number of latent variables to use in PLS classifier    
ParamsDict['NumLatents'] = 3

# Specify number of shuffles to perform to generate model test performance 
# bootstrap distributions
ParamsDict['NumRepetitions'] = 30

# Right arm to both right and left targets
ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set peri-event extraction window (seconds)
ParamsDict['BoundaryWindow'] = [-2., 2.]

# Set width of sliding window (seconds)
ParamsDict['WindowWidth'] = 0.4

# Set increment size to slide window
ParamsDict['StepWidth'] = 0.1

# Set confidence interval
ParamsDict['ConfLevel'] = 0.95

# Run batch processor
DecodingStabilityOutput = \
    DecodingStabilityBatchProcessor(ParamsDict, 
                                    long_reg_sessions_list=LongSessionsList,
                                    decoding_sessions_list=DecodingSessionsList,
                                    analysis_type='sliding_window',
                                    already_parsed=True)
    
# Save workspace
# Specify the variables to save that reside in present workspace
VarsToSave = ['DecodingStabilityOutput', 
              'DecodingSessionsList', 
              'LongSessionsList', 
#              'NumLatents', 
#              'NumShuffles', 
              'ParamsDict', 
              'SavePath']


# Specify the full path to save the shelved workspace
SavePath = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/DecodingStabilityOverPeriReachTime_20181217-18_3latents_3'
#SavePath = '/media/thugwithyoyo/STORE N GO/LongRegBatchProcessing/DecSet_LeftHemRightArm_5latents_3'

# Initialize and open self object
my_shelf = shelve.open(SavePath)

# Iterate through the list of keys in the save list, writing each to shelf
# object
for key in VarsToSave:

    # Surround write command in try/except clause to bypass possible TypeError
    # issues.
    try:

#        my_shelf[key] = globals()[key]
        my_shelf[key] = locals()[key]

    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))

# Close shelf object after variables have been written to file
my_shelf.close()
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


def DecodingStabilityBatchProcessor(NumLatents, NumShuffles, ParamsDict, **kwargs):
    
    global CellFluorTraces_Frame
    global BehavDict
    
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
    PathToDataDir = '/media/thugwithyoyo/SP UFD U3/LongRegBatchProcessing/Datafiles'
    
    # Specify path to directory that will contain shelved workspaces for each
    # session after parsing.
    PathToParsedDataDir = '/media/thugwithyoyo/SP UFD U3/LongRegBatchProcessing/ParsedWorkspaces'

    # Pack into a dict the two lists that specify the names and assigned 
    # numerical values of the behavioral event types listed in ParamsDict.
    # This dictionary will be used in later calls to the PLS decoder.
    RefEventsDict = {'RefEventsList' : ParamsDict['RefEventsList'],
                     'AssignedEventVals' : ParamsDict['AssignedEventVals']}
    
    #### Involved file handling here NOT FINISHED YET ###
    # Examine directory containing datafiles for each of the sessions in list.
    # Compare with unique sessions list from above.

    
    #os.chdir(PathToDataDir)
    onlyfiles = [f for f in listdir(PathToDataDir) if isfile(join(PathToDataDir, f))]
    
    # For each entry in the unique sessions list that is NOT in the directory
    # list, prompt user to navigate to the corresponding datafile.
    
    ####################################
    
    ########### BEGIN PARSING SESSION DATA INTO SHELVED WORKSPACES ############
    # Assume all corresponding datafiles are in data directory.
    # Iterate through each entry in  SessionsList, build path to the datafiles
    # for constructing behavioral and fluorescence data structures.  Shelve
    # datastructures to workspace files to be opened and used as necessary.
    #if kwargs['already_parsed'] is (False or None):
    if (~np.isin('already_parsed', np.array(list(kwargs.keys())))) or (kwargs['already_parsed'] == False):
        
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
                    
                    #my_shelf[key] = globals()[key]
                    my_shelf[key] = locals()[key]
                    
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
    
    NumOfCommonCellsMat = np.empty((NumSessions, NumSessions), dtype=int)
    
    TrainingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)
    
    TestingSetPerfMat = np.empty((NumSessions, NumSessions), dtype=float)
    
    ShuffTrainingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                      dtype=float)
    
    ShuffTestingSetPerfMat = np.empty((NumSessions, NumSessions, NumShuffles),
                                      dtype=float)
    
    #CellFluorTraces_Frame = pd.DataFrame()
    #BehavDict = {}
    
    ############ END STORAGE ARRAY INITIALIZATION BLOCK ###################
    
    
    ####################  BEGIN BATCH PROCESSING LOOP #####################
    SessionsListIndices = np.arange(0, NumSessions)
    ShufflesListIndices = np.arange(0, NumShuffles)
    
    # Iterate over all possible combinations of training and testing sets
    for i in SessionsListIndices:
        
        # Retrieve parsed workspaces of current training and testing
        # datasets

        # Open TRAINING set parsed workspace
        RestoreFile = SessionsList[i] + '_SessionData.dat'
        RestoreFilePath = join(PathToParsedDataDir, RestoreFile)
        exec(open('./RestoreShelvedWorkspaceScript.py').read())

        CellFluorTraces_Frame_train = CellFluorTraces_Frame         
        
        # Assign behavioral data dict for TRAINING set
        BehavDict_train = BehavDict
            
        for j in SessionsListIndices:
            
            # Obtain list of common cells between current training set 
            # (index i) and testing set (index j).
            CommonCellsDict = LRFL.CommonCellsExtractor(lr_frame, 
                                        LongRegSessionsList[i], 
                                        LongRegSessionsList[j])

            # Include only cells listed in common cells dictionary that are also
            # listed in the cell traces dataframe for training.
            Mask1 = np.isin(CommonCellsDict['S1_local_cellset'], 
                            CellFluorTraces_Frame_train.columns.values)
                        
            # Open TESTING set parsed workspace
            RestoreFile = SessionsList[j] + '_SessionData.dat'
            RestoreFilePath = join(PathToParsedDataDir, RestoreFile)
            exec(open('./RestoreShelvedWorkspaceScript.py').read())
            
            # Include only cells listed in common cells dictionary that are also
            # listed in the cell traces dataframe for testing.
            Mask2 = np.isin(CommonCellsDict['S2_local_cellset'], 
                            CellFluorTraces_Frame.columns.values)

            # Set intersect the two masks to obtain the final selection mask for
            # the cell fluorescence sets for both the training and testing sessions
            SelectFilt = (Mask1 & Mask2)
            
            # Count the number of common cells to be used for this combination
            # of training and testing sessions.
            NumOfCommonCellsMat[i, j] = np.sum(SelectFilt)
            
            # Select common cell traces from TRAINING set            
            CellFluorTraces_Frame_train = CellFluorTraces_Frame_train[
                    np.hstack([np.array(['Timestamps']), 
                               CommonCellsDict['S1_local_cellset'][SelectFilt]])]
    
            # Select common cell traces from TESTING set
            CellFluorTraces_Frame_test = CellFluorTraces_Frame[
                    np.hstack([np.array(['Timestamps']), 
                               CommonCellsDict['S2_local_cellset'][SelectFilt]])]

            # Assign behavioral data dict for TESTING set
            BehavDict_test = BehavDict
            
            # Run peri-event extraction and assemble the dictionary to be 
            # used for PLS-regression training.
            PeriEventExtractorDict_train = \
                PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                               CellFluorTraces_Frame_train, 
                                               RefEventsDict, 
                                               ParamsDict['BoundaryWindow'])
                
            # Run peri-event extraction and assemble the dictionary to be 
            # used to test the PLS-regression model above.
            PeriEventExtractorDict_test = \
                PETFL.PeriEventExtractor_Trace(BehavDict_test, 
                                               CellFluorTraces_Frame_test, 
                                               RefEventsDict, 
                                               ParamsDict['BoundaryWindow'])
                
            # Compute the PLS regression matrices from the TRAINING dataset.
            Performance_train = \
                PETFL.PLS_DecoderPerformance(PeriEventExtractorDict_train, 
                                             NumLatents)
                
            TrainingSetPerfMat[i, j] = Performance_train['performance']
            
            # Apply the test dataset to the current model trained above
            TestingSetPerfMat[i, j] = \
                PETFL.ApplyPLSModel(Performance_train['B'], 
                                    Performance_train['B_0'], 
                                    PeriEventExtractorDict_test)['performance']
                
            for k in ShufflesListIndices:
                
                # Shuffle TRAINING set cell order and run model on the result
                # Shuffle order of cell fluorescence traces in the dataframe
                CellFluorTraces_Frame_train_shuf = \
                    CellFluorTraces_Frame_train[
                            np.hstack([np.array(['Timestamps']), 
                                       np.random.choice(CommonCellsDict['S1_local_cellset'][SelectFilt], 
                                                        size=NumOfCommonCellsMat[i, j], 
                                                        replace=False)])]
                
                # Generate test data structure with shuffled cell order
                # of the TRAINING set
                PeriEventExtractorDict_train_shuf = \
                    PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                    CellFluorTraces_Frame_train_shuf, 
                                    RefEventsDict, 
                                    ParamsDict['BoundaryWindow'])
    
                # Run model on shuffled data structure and record performance
                ShuffTrainingSetPerfMat[i, j, k] = \
                    PETFL.ApplyPLSModel(Performance_train['B'], 
                                        Performance_train['B_0'], 
                                        PeriEventExtractorDict_train_shuf)['performance']
                    
                # Shuffle TESTING set cell order and run model on the result
                # Shuffle order of cell fluorescence traces in the dataframe
                CellFluorTraces_Frame_test_shuf = \
                    CellFluorTraces_Frame_test[
                            np.hstack([np.array(['Timestamps']), 
                                       np.random.choice(CommonCellsDict['S2_local_cellset'][SelectFilt], 
                                                        size=NumOfCommonCellsMat[i, j], 
                                                        replace=False)])]
                
                # Generate test data structure with shuffled cell order
                # of the TESTING set
                PeriEventExtractorDict_test_shuf = \
                    PETFL.PeriEventExtractor_Trace(BehavDict_train, 
                                    CellFluorTraces_Frame_test_shuf, 
                                    RefEventsDict, 
                                    ParamsDict['BoundaryWindow'])
    
                # Run model on shuffled data structure and record performance
                ShuffTestingSetPerfMat[i, j, k] = \
                    PETFL.ApplyPLSModel(Performance_train['B'], 
                                        Performance_train['B_0'], 
                                        PeriEventExtractorDict_test_shuf)['performance']
                        
    return {
            'NumOfCommonsCellsMat': NumOfCommonCellsMat,
            'TrainingSetPerfMat': TrainingSetPerfMat,
            'TestingSetPerfMat': TestingSetPerfMat,
            'ShuffTrainingSetPerfMat': ShuffTrainingSetPerfMat,
            'ShuffTestingSetPerfMat': ShuffTestingSetPerfMat
            }
    
    ############################# END FUNCTION ################################
    
############################# START EXECUTION SCRIPT ##########################
    
LongSessionsList = np.array(['2018-12-04-10-31-21',
                             '2018-12-10-11-37-56',
                             '2018-12-11-10-53-04',
                             '2018-12-14-11-01-41',
                             '2018-12-17-11-38-42',
                             '2018-12-18-11-20-21',
                             '2018-12-28-10-50-15'])
    
DecodingSessionsList = np.array(['2018-12-04-10-31-21',
                                 '2018-12-10-11-37-56',
                                 '2018-12-11-10-53-04',
                                 '2018-12-14-11-01-41',
                                 '2018-12-17-11-38-42',
                                 '2018-12-18-11-20-21',
                                 '2018-12-28-11-12-14'])

# Specify number of latent variables to use in PLS classifier    
NumLatents = 5

# Specify number of shuffles to perform to generate model test performance 
# bootstrap distributions
NumShuffles = 2
    
# Define a dictionary to contain keyed access to sliding window analysis parameters
ParamsDict = defaultdict(dict)

# Right arm to both right and left targets
ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']

# Scalar values assigned to event types listed above.
ParamsDict['AssignedEventVals'] = [-1, 1]

# Set peri-event extraction window (seconds)
ParamsDict['BoundaryWindow'] = [-1., 1.]

# Run batch processor
DecodingStabilityDict = \
    DecodingStabilityBatchProcessor(NumLatents, NumShuffles, ParamsDict, 
                                    long_reg_sessions_list=LongSessionsList,
                                    decoding_sessions_list=DecodingSessionsList,
                                    already_parsed=True)
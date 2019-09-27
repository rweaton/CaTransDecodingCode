#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 09:19:42 2019

@author: thugwithyoyo
"""
import numpy as np
import pandas as pd
from collections import defaultdict

# Example calling line.
#    OutputTraceMatrix, ProcessingDict = CaTraceNormalizer(
#                                              AveragedTraceMatrices[:,:,i], 
#                                              CellFluorTraces_Frame,
#                                              ParamsDict,
#                                              NormalizationMethod)

def zScoreTraces(FullTraces, ParamsDict):
    
    # Talk about condescending..!  You're research is a bit of a joke.  Hers is
    # for real.  
    (NumTraces, NumSamples) = FullTraces.shape
    
    TraceMeans = np.array([np.mean(FullTraces, axis=1)]).transpose()
    TraceStdDevs = np.array([np.std(FullTraces, axis=1)]).transpose()
    
    RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
                             ParamsDict['BoundaryWindow'][1], 
                             num=NumSamples, endpoint=False)
    
    # Z-score the rows of MatrixOfTraces using corresponding means and standard
    # deviations computed from the complete records.
    
    zScoredTraces = np.divide(
        (FullTraces - np.dot(TraceMeans,  np.ones((1, NumSamples), dtype=float))),
        np.dot(TraceStdDevs, np.ones((1, NumSamples), dtype=float)))

    return zScoredTraces

def CaTraceSorter(MatrixOfTraces, ParamsDict, SortMethod):
    
    # MatrixOfTraces must contain each trace as a row.
    
    # Specify the domain, relative to restrict search for peaks or other features.
    #SearchDomain = [-1., 0.]
    
    # Initialize a dictionary to contain normalization and sorting processing output.
    ProcessingDict = defaultdict(dict)
    
    # Acquire mean and standard deviation statistics from entire record from
    # which peri-reach traces were originally extracted.  This is not valid
    # since the population of traces are AVERAGES of calcium activity!!!
    #nCols = CellFluorTraces_Frame.shape[1]
    
    # Generate a row matrix array of full-length ca trace records.
    #FullTraces = pd.DataFrame(CellFluorTraces_Frame._slice(slice(1, nCols), 1)).values.transpose()

    # Generate a row vector array of the list of time corresponding to 
    # full-length records.
    #Timestamps_df = pd.DataFrame(CellFluorTraces_Frame._slice(slice(0, 1), 1)).values.transpose()
    

        
    # Column vectors of means and standard deviations of each trace.
    #TraceMeans = np.array([np.mean(FullTraces, axis=1)]).transpose()
    #TraceStdDevs = np.array([np.std(FullTraces, axis=1)]).transpose()
    #TraceMeans = np.array([np.mean(MatrixOfTraces, axis=1)]).transpose()
    #TraceStdDevs = np.array([np.std(MatrixOfTraces, axis=1)]).transpose()
    
    (NumTraces, NumSamples) = MatrixOfTraces.shape
    RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
                             ParamsDict['BoundaryWindow'][1], 
                             num=NumSamples, endpoint=False)
    
    # Z-score the rows of MatrixOfTraces using corresponding means and standard
    # deviations computed from the complete records.
#        zScoredTraces = np.divide(
#            (MatrixOfTraces - np.dot(TraceMeans,  np.ones((1, NumSamples), dtype=float))),
#            np.dot(TraceStdDevs, np.ones((1, NumSamples), dtype=float)))
    
    # Setup a filter for defining scope to identify feature search (e.g. peak level) 
    SearchDomainFilt = (RelTimeVec >= ParamsDict['SearchDomain'][0]) & \
                       (RelTimeVec <= ParamsDict['SearchDomain'][1])
    
    SearchDomainIndices = np.tile(np.array([np.arange(0, SearchDomainFilt.shape[0], 1)]).transpose(), 
                                    NumTraces).transpose()
    if SortMethod == 'PeakSort':
        
        # Locate maxima within search domain and computer their values.    
        ProcessingDict['MaxLocs'] = np.array([np.argmax(MatrixOfTraces[:, SearchDomainFilt], 
                                             axis=1)]).transpose()
    
    elif SortMethod == 'AbsPeakSort':
        
        # Locate maxima within search domain and computer their values.    
        ProcessingDict['MaxLocs'] = np.array([np.argmax(np.abs(MatrixOfTraces[:, SearchDomainFilt]), 
                                             axis=1)]).transpose()        
    
    ValExtractionFilt = (SearchDomainIndices == ProcessingDict['MaxLocs']) 
    
    ProcessingDict['MaxVals'] = MatrixOfTraces[ValExtractionFilt]
    
    # Generate a sorting map to organize row traces.
    ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['MaxVals'], axis=-1)
    
    # Sort traces       
    #OutputTraceMatrix = MatrixOfTraces[ProcessingDict['SortIndices'], :]
           
        
    #return OutputTraceMatrix, ProcessingDict
    return ProcessingDict

def PeakNormalizer(MatrixOfTraces):
            
    (NumTraces, NumSamples) = MatrixOfTraces.shape
    TraceAbsValMaxima = np.max(np.abs(MatrixOfTraces), axis=-1)
    
    NormVecs = np.dot(np.array([TraceAbsValMaxima]).transpose(), np.array([np.ones((NumSamples,))]))
    
    return np.divide(MatrixOfTraces, NormVecs) 
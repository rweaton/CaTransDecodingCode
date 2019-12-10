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
        
        # Locate trace maxima within search domain and computer their values.    
        ProcessingDict['MaxLocs'] = np.array([np.argmax(MatrixOfTraces[:, SearchDomainFilt], 
                                             axis=1)]).transpose()
    
        ValExtractionFilt = (SearchDomainIndices == ProcessingDict['MaxLocs']) 
    
        ProcessingDict['MaxVals'] = MatrixOfTraces[ValExtractionFilt]
        
        ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['MaxVals'], axis=-1)
    
    elif SortMethod == 'AbsPeakSort':
        
        # Locate maxima of absolute valued traces within search domain and computer their values.    
        ProcessingDict['MaxLocs'] = np.array([np.argmax(np.abs(MatrixOfTraces[:, SearchDomainFilt]), 
                                             axis=1)]).transpose()
        ValExtractionFilt = (SearchDomainIndices == ProcessingDict['MaxLocs']) 
    
        ProcessingDict['MaxVals'] = MatrixOfTraces[ValExtractionFilt]
        
        ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['MaxVals'], axis=-1)
    
    elif SortMethod =='AreaSort':

        ProcessingDict['AreaVals'] = np.array([np.sum(MatrixOfTraces[:, SearchDomainFilt],
                                              axis=1)]).transpose()
        
        ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['AreaVals'], axis=0).transpose()[0]
        
    elif SortMethod == 'AvgSort':
        
        ProcessingDict['AvgVals'] = np.array([np.mean(MatrixOfTraces[:, SearchDomainFilt],
                                              axis=-1)]).transpose()
    
        ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['AvgVals'], axis=0).transpose()[0]
 
#    ValExtractionFilt = (SearchDomainIndices == ProcessingDict['MaxLocs']) 
    
#    ProcessingDict['MaxVals'] = MatrixOfTraces[ValExtractionFilt]
    
    # Generate a sorting map to organize row traces.
#    ProcessingDict['SortIndices'] = np.argsort(ProcessingDict['MaxVals'], axis=-1)
    
    # Sort traces       
    #OutputTraceMatrix = MatrixOfTraces[ProcessingDict['SortIndices'], :]
           
        
    #return OutputTraceMatrix, ProcessingDict
    return ProcessingDict

def GetMax(Array):
    
    if len(Array.shape) > 1:
        
        return GetMax(np.amax(Array, axis=-1))

    else:
        
        return np.amax(Array)

def GetMean(Array):

    if len(Array.shape) > 1:
        
        return GetMean(np.mean(Array, axis=-1))

    else:
        
        return np.mean(Array)
    
def GetTuning(ArrayOfMatrices, TuningMethod, **kwargs):
    
    DifferenceMatrix = np.diff(ArrayOfMatrices, axis=-1)

    if TuningMethod == 'MaxAbsDiffMatrixNorm':
        
        NormVal = GetMax(np.abs(DifferenceMatrix))
        
        TuningIndexArray = DifferenceMatrix/NormVal
        
    if TuningMethod == 'DiffOverSum':
#    NormVal = GetMean(np.abs(DifferenceMatrix))
        
        TraceAvgs = np.mean(ArrayOfMatrices, axis=1)
        
        NormVal = np.sum(np.abs(TraceAvgs), axis=-1)
        
        TuningIndexArray = np.divide(TraceAvgs[:,1] - TraceAvgs[:,0],
                                     NormVal)
        
    if TuningMethod == 'DiffOfAvgsOverSumOfMagOfAvgs':
        
        # Begin calculation of scalar tuning indices for z-score averaged traces from 
        # each cell
        (NumTraces, NumSamples) = ArrayOfMatrices[:,:,0].shape
        
        RelTimeVec = np.linspace(kwargs['ParamsDict']['BoundaryWindow'][0], 
                                 kwargs['ParamsDict']['BoundaryWindow'][1], 
                                 num=NumSamples, endpoint=False)
        
        SearchDomainFilt = (RelTimeVec >= kwargs['ParamsDict']['SearchDomain'][0]) & \
                           (RelTimeVec <= kwargs['ParamsDict']['SearchDomain'][1])
                           
        TraceAvgs = np.mean(ArrayOfMatrices[:, SearchDomainFilt, :], axis=1)
        
        NormVal = np.sum(np.abs(TraceAvgs), axis=-1)
        
        TuningIndexArray = np.divide(TraceAvgs[:,1] - TraceAvgs[:,0],
                                     NormVal)
        
    if TuningMethod == 'WeightedDifference':
        
        # Begin calculation of scalar tuning indices for z-score averaged traces from 
        # each cell
        (NumTraces, NumSamples) = ArrayOfMatrices[:,:,0].shape
        
        RelTimeVec = np.linspace(kwargs['ParamsDict']['BoundaryWindow'][0], 
                                 kwargs['ParamsDict']['BoundaryWindow'][1], 
                                 num=NumSamples, endpoint=False)
        
        SearchDomainFilt = (RelTimeVec >= kwargs['ParamsDict']['SearchDomain'][0]) & \
                           (RelTimeVec <= kwargs['ParamsDict']['SearchDomain'][1])
                           
        TraceAvgs = np.mean(ArrayOfMatrices[:, SearchDomainFilt, :], axis=1)
        
        #mu = GetMean(ArrayOfMatrices)*np.ones(TraceAvgs.shape[0])
        
        Weights = np.exp(-(1./2.)*((np.sum(TraceAvgs, axis=-1))/kwargs['ParamsDict']['StdDev'])**2)
        
        TuningIndexArray = np.multiply((TraceAvgs[:,1] - TraceAvgs[:,0]), Weights)
        
    if TuningMethod == 'PlainDifference':
        
        # Begin calculation of scalar tuning indices for z-score averaged traces from 
        # each cell
        (NumTraces, NumSamples) = ArrayOfMatrices[:,:,0].shape
        
        RelTimeVec = np.linspace(kwargs['ParamsDict']['BoundaryWindow'][0], 
                                 kwargs['ParamsDict']['BoundaryWindow'][1], 
                                 num=NumSamples, endpoint=False)
        
        SearchDomainFilt = (RelTimeVec >= kwargs['ParamsDict']['SearchDomain'][0]) & \
                           (RelTimeVec <= kwargs['ParamsDict']['SearchDomain'][1])
                           
        TraceAvgs = np.mean(ArrayOfMatrices[:, SearchDomainFilt, :], axis=1)
        
        TuningIndexArray = (TraceAvgs[:,1] - TraceAvgs[:,0])
        
    return TuningIndexArray


def PeakNormalizer(MatrixOfTraces):
            
    (NumTraces, NumSamples) = MatrixOfTraces.shape
    TraceAbsValMaxima = np.max(np.abs(MatrixOfTraces), axis=-1)
    
    NormVecs = np.dot(np.array([TraceAbsValMaxima]).transpose(), np.array([np.ones((NumSamples,))]))
    
    return np.divide(MatrixOfTraces, NormVecs) 

def TuningByCellFrameGen(CellFluorTraces_Frame, SortProcessingDict, ParamsDict):
    
    # Write into dataframe: cell identities, tuning values and indicies 
    # of traces as they are tuning ordered in the heatmaps.
    
    # Count the total number of columns in the fluorescence trace dataframe
    (_, NumColumns) = CellFluorTraces_Frame.shape

    # Number of cells included is one less than the total number of columns.  
    # The first column contains timestamps.
    NumCells = NumColumns - 1
    
    # Extract the cell names from the column header.
    CellLabels = np.array(list(CellFluorTraces_Frame.columns.values))[1:]
    
    # Initialize tuning dataframe.
    Tuning_Frame = pd.DataFrame(index=CellLabels, 
                                columns=['ScalarTuningIndex', 'HeatMapRowIndex'])
    
    # Write tuning indices to the tuning dataframe.
    Tuning_Frame['ScalarTuningIndex'] = SortProcessingDict['AvgVals'].transpose()[0]
    
    # Perform reverse mapping operation to indicate each the row index of the 
    # associated trace in the tuning-sorted heat map.
    
    # Generate a list of indices for the cells
    IndexList = np.arange(0,NumCells)
    
    # Iterate through the list of indices and detect location of same index in
    # the list of indices sorted by corresponding tuning index value.
    for i in IndexList:
        
        # Find location of trace index in the sorted list.
        Filt = (SortProcessingDict['SortIndices'] == i) 
        
        # Write CellLabel and index location to the tuning dataframe.
        Tuning_Frame.at[CellLabels[i], 'HeatMapRowIndex'] = IndexList[Filt][0]
    
    return Tuning_Frame
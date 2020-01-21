#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 05:53:13 2019

@author: thugwithyoyo
"""

import json
import os
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Access examples
#TiltDataDict['neurons']['sig003a']
#TiltDataDict['neurons']['sig016b']
#TiltDataDict['events']['event_6']
#TiltDataDict['events']['event_3']


def PeriEventTimeStampExtractor(EventTimesList, TimeStampsList, BoundaryWindow):
    
    EventTimesList = np.array(EventTimesList)
    TimeStampsList = np.array(TimeStampsList)
    
    eIndices = np.arange(0, EventTimesList.size)
    PeriEventTimes = EventTimesList.size*[None]
    
    for i in eIndices:
        # This is only one of four possible extraction definitions.  Will use
        # this as the default.  Consider including a case construct with a 
        # an optional flag if the user would like to use one of the other
        # possible definitions.
        filt = ((TimeStampsList >= EventTimesList[i] + BoundaryWindow[0]) &  
                (TimeStampsList < EventTimesList[i] + BoundaryWindow[1]))
        
#        filt = ((TimeStampsList > EventTimesList[i] + BoundaryWindow[0]) &  
#                (TimeStampsList <= EventTimesList[i] + BoundaryWindow[1]))        
       
        PeriEventTimes[i] = np.array(TimeStampsList[filt] - EventTimesList[i])
        
    return PeriEventTimes

def RasterMatrixCompiler(PeriEventTimes, BoundaryWindow, Resolution):
    
    # Argument must be a list of arrays of event relative list of times to be
    # compiled into a binary matrix; each row indicating trial no. and each
    # each column entry in that row indicating the presence (1) or absence (0)
    # of an event within a bin of width Resolution.
    
    # Number of Events
    nEvents = len(PeriEventTimes)
    eIndices = np.arange(0, nEvents)
    
    # Max spike resolution
    #Resolution = 0.0005
    
    # Pre-allocate the binary occurance matrix.
    # Compute number of boxes needed from Resolution
    nBins = int((BoundaryWindow[1] - BoundaryWindow[0])/Resolution)
    
    # The occurrence matrix is size [nEvents X nBins]
    RasterMatrix = np.nan*np.ones([nEvents, nBins])
    RasterMatrix_b = np.zeros([nEvents, nBins])
    #TimeBinVec = np.linspace(BoundaryWindow[0], BoundaryWindow[1], 
    #                         num=nBins, endpoint=True)
    #TimeBinVec = np.arange(BoundaryWindow[0]/Resolution,  
    #                       (BoundaryWindow[1]/Resolution)
    BinBoundaryVec = np.linspace(BoundaryWindow[0]/Resolution, 
                             BoundaryWindow[1]/Resolution, 
                             num = (nBins+1), 
                             endpoint = True)
    
    BinBoundaryVec = np.around(BinBoundaryVec)
    
    BinCenters = ((BinBoundaryVec[1:] - 
                   BinBoundaryVec[0:-1])/2. + BinBoundaryVec[0:-1])
    
    nFiltEvents = np.nan*np.ones_like(eIndices)
    for i in eIndices:
        
        filt = np.isin(BinBoundaryVec[0:-1], np.floor(PeriEventTimes[i]/Resolution), 
                       assume_unique=False)
        
        nFiltEvents[i] = sum(filt)
        
        RasterMatrix[i, filt] = Resolution*BinCenters[filt]
        RasterMatrix_b[i, filt] = np.ones_like(RasterMatrix_b[i,:])[filt]
    
    return {
            'Resolution' : Resolution,
            'nFiltEvents' : nFiltEvents,
            'RasterMatrix' : RasterMatrix,
            'RasterMatrix_b' : RasterMatrix_b,
            'BinBoundaryVec' : Resolution*BinBoundaryVec,
            'BinLeftEdgesVec' : Resolution*BinBoundaryVec[0:-1],
            'BinCentersVec' : Resolution*BinCenters
           }

def PeriEventHistoFromRaster(RasterMatrix, BinBoundaryVec):
    # NOTE: This routine is only valid for bin-widths that capture only
    # a single event per bin!!!  Use PeriEventHistoFromRaster2 for wider bins.
    (nTrials, nBins) = RasterMatrix.shape
    OccurrenceMatrix = ~np.isnan(RasterMatrix)
    Counts = np.sum(OccurrenceMatrix, axis=0)
    
    BinCenters = ((BinBoundaryVec[1:] - 
                       BinBoundaryVec[0:-1])/2. + BinBoundaryVec[0:-1])
    
    BinWidth = BinBoundaryVec[1] - BinBoundaryVec[0]
    
    return {
            'nTrials':nTrials,
            'Counts':Counts,
            'NormalizedHistogram' : Counts/nTrials,
            'BinBoundaryVec': BinBoundaryVec,
            'BinLeftEdgesVec' : BinBoundaryVec[0:-1],
            'BinCentersVec' : BinCenters
            }

def PeriEventHistoFromRaster2(RasterDict, HistoBinWidth):
    
    (nTrials, nRasterBins) = RasterDict['RasterMatrix'].shape
    
    # HistoBinWidth should be an integer multiple of raster resolution
    nRasterBinsPerHistoBin = int(np.floor(HistoBinWidth/RasterDict['Resolution'] )) # Uncertain if floor operation is valid
    
    SlicingArray = np.hstack((np.arange(0, nRasterBins, nRasterBinsPerHistoBin), 
                             np.array([nRasterBins])))
    
    Counts = np.nan*np.ones_like(SlicingArray[0:-1])
    
    for i in np.arange(0, (SlicingArray.size - 1)):
        
        ColumnSliceDef = slice(SlicingArray[i], SlicingArray[i + 1], 1)
        Counts[i] = np.sum(np.sum(RasterDict['RasterMatrix_b'][:, ColumnSliceDef], axis=1), axis=0)
    
    #OccurrenceMatrix = ~np.isnan(RasterMatrix)
    #Counts = np.sum(OccurrenceMatrix, axis=0)
    
    BinBoundaryVec = RasterDict['BinBoundaryVec'][SlicingArray]
    
    BinCenters = ((BinBoundaryVec[1:] - 
                   BinBoundaryVec[0:-1])/2. + BinBoundaryVec[0:-1])
    
    return {
            'nTrials':nTrials,
            'Counts':Counts,
            'NormalizedHistogram' : Counts/nTrials,
            'BinBoundaryVec': BinBoundaryVec,
            'BinLeftEdgesVec' : BinBoundaryVec[0:-1],
            'BinCentersVec' : BinCenters
            }
    
def ReceptiveFieldAnalysis(PETH_Dict, Thresh_nSD, nContigBins):
    
    BinIndicesVec = np.arange(0, PETH_Dict['BinLeftEdgesVec'].size)
    #BinWidth = PETH_Dict['BinLeftEdgesVec'][1] - PETH_Dict['BinLeftEdgesVec'][0]
    
    ### Compute stats on spontaneous activity
    # Define domain of spont. activity
    filt = PETH_Dict['BinLeftEdgesVec'] < 0.
    
    # Calculate threshold for post event response.
    background_rate = np.mean(PETH_Dict['NormalizedHistogram'][filt])
    background_SD = np.std(PETH_Dict['NormalizedHistogram'][filt], ddof=1)
    threshold = background_rate + Thresh_nSD*background_SD
    
    ### Locate and measure post-event effect
    # Extract post-event domain
    filt1 = PETH_Dict['BinLeftEdgesVec'] >= 0.
    
    # Extract superthreshold bins
    filt2 = (PETH_Dict['NormalizedHistogram'] > threshold)
    
    # Slide over and cumulatively sum traces (slid forward and slid backward)  
    # of above-theshold binary filter traces.
    ContigBinTrace_fw = GenContigBinTrace(filt2, nContigBins, 'forward')
    ContigBinTrace_bw = GenContigBinTrace(filt2, nContigBins, 'backward')
    
    # Use the contiguous bins trace to find the post event effect.
    # Find maximum of post event effect 
    MaxContigBins = np.max(ContigBinTrace_fw)
    EffectMaximaFilt = ((ContigBinTrace_fw == MaxContigBins) & filt1)
    
    # For each maxima in the set, trace left and right to find respective
    # effect start and effect end locations. Add values to lists then filter
    # them for uniqueness.
    EffectMaximaIndices = BinIndicesVec[EffectMaximaFilt]
    EffectStartIndices = []
    EffectEndIndices = []
    
    for i in EffectMaximaIndices:

        # Trace left to first zero encountered.  This is when the effect
        # first fails the contiguity requirement; indicating effect onset.        
        j = i
        
        while ((ContigBinTrace_bw[j] > 0) and (j > 0)):
            
            j -= 1
            
        #EffectStartIndices.append(j + 1)
        EffectStartIndices.append(j)

        # Trace right to first zero encountered.  This is when the effect
        # first fails the contiguity requirement to indicate effect 
        # termination.  The effect may stay above threshold past the end of
        # the window so the trace loop must account for this end case.        
        j = i
    
        while ((ContigBinTrace_fw[j] > 0) and (j <= ContigBinTrace_fw.size - 2)):
            
            j += 1
        
        # Check for end bin separately (nasty hack solution)
        if (PETH_Dict['NormalizedHistogram'][-1] > threshold):
        
            j += 1
            
        EffectEndIndices.append(j - 1)    
            
    EffectIndicesArray = np.array([EffectStartIndices, EffectEndIndices]).transpose()
    
    # Should remove redundant entries in EffectIndicesArray here...  Not sure 
    # of an elegant way to do it right now.  Simply going to take the first
    # element here.
    EffectIndices = EffectIndicesArray[0]
    
    #EffectIndices = BinIndicesVec[(ContigBinTrace > 0) & filt1]
    
    # Find relative time domain of the post event effect.
    first_bin_latency = PETH_Dict['BinCentersVec'][EffectIndices[0]]
    last_bin_latency = PETH_Dict['BinCentersVec'][EffectIndices[1]]
    
    #EffectDomain = np.linspace(EffectIndices[0], EffectIndices[1], endpoint=True)
    EffectDomain = slice(EffectIndices[0], EffectIndices[1] + 1, 1)
    
    # Find value and latency of response peak
    PeakVal = np.max(PETH_Dict['NormalizedHistogram'][EffectDomain])
    filt = (PETH_Dict['NormalizedHistogram'] == PeakVal)
    peak_latency = PETH_Dict['BinCentersVec'][filt]
    peak_firing_rate = PeakVal

    # Determine response magnitude
    #EffectDomain = slice(EffectIndices[0], EffectIndices[1] + 1, 1)
    response_magnitude = np.sum(PETH_Dict['Counts'][EffectDomain])/PETH_Dict['nTrials']
    psth = PETH_Dict['NormalizedHistogram']
    #filt2 = (ContigBinTrace > nContigBins) & (filt1)
    
    ReceptiveFieldDict = { 
                           'psth' : psth.tolist(),
                           'background_rate' : background_rate,
                           'background_SD' : background_SD,
                           'threshold' : threshold,
                           'ContigBinTrace_fw' : ContigBinTrace_fw.tolist(),
                           'ContigBinTrace_bw' : ContigBinTrace_bw.tolist(),
                           'first_bin_latency' : first_bin_latency,
                           'last_bin_latency' : last_bin_latency,
                           'peak_latency' : peak_latency[0],
                           'peak_firing_rate' : peak_firing_rate,
                           'response_magnitude' : response_magnitude,
                         }
    
    return ReceptiveFieldDict

def GenContigBinTrace(bSequence, nContigBins, SlideDir):
    # Function to calculate the contiguous bin count for computing effect
    # boundaries in peri-stimulus.  It works by sliding a binary "above-
    # threshold" vector over by 1 to nContigBins then summing all of the
    # resulting displaced traces together, element-by-element.  Where the
    # summed trace first falls to zero indicates effect boundaries. The binary
    # trace can be slid in either direction as specified by the SlideDir flag.
    
    # Convert "true" or "false" to "ones" or "zeros"
    bSequence = bSequence.astype(int)
    
    # Count the number of elements in the binary vector
    (nSeqBins,) = bSequence.shape
    
    # Case in which traces are sequentially slid in positive direction
    if SlideDir == 'forward':
        
        # Generate the list of integer slide displacements
        IterationList = np.arange(0, nContigBins)
        
        # Pad ends of sequence array for subsequent shift-and-add operation.
        bSequence_padded = np.pad(bSequence, (0, nContigBins), 'constant')
    
        # Initialize the result trace
        ContigBinTrace = np.zeros_like(bSequence)
        
        # Extract out a slice of the zero-padded binary trace that is displaced 
        # positively by i elements. Add the trace (element-wise) to cumulative
        # sum.
        for i in IterationList:
        
            ContigBinTrace = ContigBinTrace + bSequence_padded[i:(nSeqBins + i)]
    
    # Case in which traces are sequentially slid in negative direction
    if SlideDir == 'backward':
        
        # Generate the list of integer slide displacements
        IterationList = np.arange(0, nContigBins)
        
        # Pad ends of sequence array for subsequent shift-and-add operation.
        bSequence_padded = np.pad(bSequence, (nContigBins, 0), 'constant')
        
        # Initialize the result trace
        ContigBinTrace = np.zeros_like(bSequence)
    
        # Extract out a slice of the zero-padded binary trace that is displaced 
        # negatively by i elements. Add the trace (element-wise) to cumulative
        # sum.
        for i in IterationList:
        
            ContigBinTrace = ContigBinTrace + bSequence_padded[
                    (nContigBins - 1 - i):(bSequence_padded.size - 1 - i)]
            
    return ContigBinTrace

def RasterPlot(RastDict, axs, show_xAxis, show_yAxis):
    
    # Extract window boundaries from bin boundary vector in raster dictionary
    BoundaryWindow = (RastDict['BinBoundaryVec'][0],
                      RastDict['BinBoundaryVec'][-1])
    
    plt.sca(axs)   
    
    # Write title of axes to figure
    plt.title(e + ', ' + n)
    
    # Count number of trials to set limits of y-axis
    yLimit = RastDict['RasterMatrix'].shape[0]
    
    # Plot raster on specified axes
    axs.eventplot(RastDict['RasterMatrix'], color='black')

    # Plot event reference line at zero
    axs.plot([0,0], [0, yLimit], linewidth=0.5, color='gray')
    
    # Hide x and y axes unless flags indicate to  show them
    plt.xlim(BoundaryWindow)
    axs.get_xaxis().set_visible(False)
    plt.ylim((yLimit, 0))
    axs.get_yaxis().set_visible(False)
    
    if (show_xAxis == True):
        axs.get_xaxis().set_visible(True)
        plt.xlabel('time (s)')
    
    if (show_yAxis == True):
        axs.get_yaxis().set_visible(True)
        plt.ylabel('trial no.')
    
def HistogramPlot(RastDict, PEHistoDict, RF_Dict, axs, show_xAxis, show_yAxis,
                  show_RFA):
    
    # Extract window boundaries from bin boundary vector in raster dictionary
    BoundaryWindow = (RastDict['BinBoundaryVec'][0],
                      RastDict['BinBoundaryVec'][-1])
    
    # Determine bin width from first two entries from bin edges vector; assumes
    # uniform width  among all bins
    BinWidth = PEHistoDict['BinLeftEdgesVec'][1] - PEHistoDict['BinLeftEdgesVec'][0]    
    
    # Plot histogram bars
    yLimit = np.ceil(20*1.05*np.max(PEHistoDict['NormalizedHistogram']))/20.
    plt.sca(axs)
    
    # Write axes title using event and neuron source identities 
    plt.title(e + ', ' + n)
    
    # Use bar graph to plot normalized histogram
    axs.bar(PEHistoDict['BinLeftEdgesVec'], 
                PEHistoDict['NormalizedHistogram'], 
                BinWidth, align='edge', color='black')
    
    # Plot reference line (zero-vertical)
    axs.plot([0, 0], [0, yLimit], linewidth=0.5, color='gray')
    
    # If show_RFA set to true, plot threshold, first bin and last bin lines
    if show_RFA:
    
        # Plot threshold line
        axs.plot(BoundaryWindow, 
                     [RF_Dict['threshold'], RF_Dict['threshold']],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
                    
        # Plot onset of post-event effect
        axs.plot([RF_Dict['first_bin_latency'], RF_Dict['first_bin_latency']], 
                     [0, yLimit],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
        
        # Plot end of post-event effect
        axs.plot([RF_Dict['last_bin_latency'], RF_Dict['last_bin_latency']], 
                     [0, yLimit],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
    
    # Hide axes numbers and labels unless specified by show_?Axis flags
    plt.xlim(BoundaryWindow)
    axs.get_xaxis().set_visible(False)
    plt.ylim((0, yLimit))
    axs.get_yaxis().set_visible(False)
    
    if (show_xAxis == True):
        axs.get_xaxis().set_visible(True)
        plt.xlabel('time (s)')
    
    if (show_yAxis == True):
        axs.get_yaxis().set_visible(True)
        plt.ylabel('spikes/bin per stim.')    
    
def RasterOverHistogramPlot(RastDict, PEHistoDict, RF_Dict, show_RFA):
    
    # Extract window boundaries from bin boundary vector in raster dictionary
    BoundaryWindow = (RastDict['BinBoundaryVec'][0],
                      RastDict['BinBoundaryVec'][-1])
    
    # Determine bin width from first two entries from bin edges vector; assumes
    # uniform width  among all bins
    BinWidth = PEHistoDict['BinLeftEdgesVec'][1] - PEHistoDict['BinLeftEdgesVec'][0]
    
    # Plot raster on upper subplot
    fig, axs = plt.subplots(nrows=2, ncols=1)
    
    # select top axes plane
    plt.sca(axs[0])
    
    # Write title of subplot
    #plt.title(n + ' spike activity aligned by ' + e)
    plt.title(e + ', ' + n + ' (bin-width = ' + str(np.round(BinWidth, decimals=3)) + ' sec.)')
    yLimit = RastDict['RasterMatrix'].shape[0]
    
    # Plot raster with axes labels in upper subplot
    axs[0].eventplot(RastDict['RasterMatrix'], color='black')
    axs[0].plot([0,0], [0, yLimit], linewidth=0.5, color='gray')
    plt.sca(axs[0])
    plt.xlim(BoundaryWindow)
    plt.ylim((yLimit, 0))
    plt.ylabel('trial no.')
    axs[0].get_xaxis().set_visible(False)
    
    # Plot histogram in lower subplot
    yLimit = np.ceil(20*1.05*np.max(PEHistoDict['NormalizedHistogram']))/20.
    plt.sca(axs[1])
    axs[1].bar(PEHistoDict['BinLeftEdgesVec'], 
                PEHistoDict['NormalizedHistogram'], 
                BinWidth, align='edge', color='black')
    
    # Plot reference line (zero-vertical)
    axs[1].plot([0, 0], [0, yLimit], linewidth=0.5, color='gray')
    
    if show_RFA:
    
        # Plot threshold line
        axs[1].plot(BoundaryWindow, 
                     [RF_Dict['threshold'], RF_Dict['threshold']],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
                    
        # Plot onset of post-event effect
        axs[1].plot([RF_Dict['first_bin_latency'], RF_Dict['first_bin_latency']], 
                     [0, yLimit],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
        
        # Plot end of post-event effect
        axs[1].plot([RF_Dict['last_bin_latency'], RF_Dict['last_bin_latency']], 
                     [0, yLimit],
                     linewidth=0.5,
                     linestyle='dashed',
                     color='red')
    
    plt.xlim(BoundaryWindow)
    plt.xlabel('time (s)')
    plt.ylim((0, yLimit))
    plt.ylabel('spikes/bin per stim.')
    
#### Begin homework 2 function definitions. ####
    
def RasterBinWidener(RasterDict, TargetBinWidth):
    
    (nTrials, nRasterBins) = RasterDict['RasterMatrix'].shape
    
    # TargetBinWidth should be an integer multiple of raster resolution
    nRasterBinsPerHistoBin = int(np.floor(TargetBinWidth/RasterDict['Resolution'] )) # Uncertain if floor operation is valid
    
    SlicingArray = np.hstack((np.arange(0, nRasterBins, nRasterBinsPerHistoBin), 
                             np.array([nRasterBins])))

    BinBoundaryVec = RasterDict['BinBoundaryVec'][SlicingArray]
    
    BinCenters = ((BinBoundaryVec[1:] - 
                   BinBoundaryVec[0:-1])/2. + BinBoundaryVec[0:-1])    
#    Counts = np.nan*np.ones_like(SlicingArray[0:-1])
    
    NewRasterMatrix_b = np.nan*np.ones([nTrials, SlicingArray[0:-1].size])
    NewRasterMatrix = np.nan*np.ones([nTrials, SlicingArray[0:-1].size])
    
    for i in np.arange(0, nTrials):
        
        for j in np.arange(0, (SlicingArray.size - 1)):
            
            ColumnSliceDef = slice(SlicingArray[j], SlicingArray[j + 1], 1)
            NewRasterMatrix_b[i, j] = np.sum(RasterDict['RasterMatrix_b'][i, ColumnSliceDef], axis=0)
    
        filt = (NewRasterMatrix_b[i, :] >= 1)
        NewRasterMatrix[i, filt] = BinCenters[filt]        

    return {
            'Resolution' : RasterDict['Resolution'],
            'TargetBinWidth' : TargetBinWidth,
#            'nFiltEvents' : nFiltEvents,
            'RasterMatrix' : NewRasterMatrix,
            'RasterMatrix_b' : NewRasterMatrix_b,
            'BinBoundaryVec' : BinBoundaryVec,
            'BinLeftEdgesVec' : BinBoundaryVec[0:-1],
            'BinCentersVec' : BinCenters
           }
    
    
def CountEntropyCalculator(RasterDict):
    
    PostStimFilt = RasterDict['BinLeftEdgesVec'] >= 0.
    
    CountsVec = np.sum(RasterDict['RasterMatrix_b'][:,PostStimFilt], axis=1)
    
    UniqueCountsSet = np.unique(CountsVec)
    
    nCountInstances = np.nan*np.ones_like(UniqueCountsSet)
    
    nCountsIndices = np.arange(0, nCountInstances.size)
    
    for i in nCountsIndices:
        
        nCountInstances[i] = np.sum(CountsVec == UniqueCountsSet[i], axis=0)
        
    CountProportions = nCountInstances/nCountInstances.sum()
    
    return {
            'Instances' : np.transpose(np.array([CountsVec])),
            'UniqueInstancesSet' : np.transpose(np.array([UniqueCountsSet])),
            'SummedInstances' : nCountInstances.sum(),
            'InstanceProportions' : CountProportions,
            'Entropy' : - np.sum(CountProportions*np.log2(CountProportions))          
            }
    
# Misunderstood what was meant by "timing entropy" here    
#def TimingEntropyCalculator_incorrect(RasterDict):
#    
#    PostStimFilt = RasterDict['BinLeftEdgesVec'] >= 0.
#    
#    BinIndices = np.arange(0, RasterDict['BinLeftEdgesVec'][PostStimFilt].size)
#    
#    OccurrenceFilt = (np.sum(RasterDict['RasterMatrix_b'][:, BinIndices], axis=0) > 0)
#    
#    OccupiedBinsSet = BinIndices[OccurrenceFilt]
#    
#    BinSetIndices = np.arange(0, OccupiedBinsSet.size)
#    
#    nBinOccupancies = np.nan*np.ones_like(OccupiedBinsSet)
#    
#    for i in BinSetIndices:
#        
#        nBinOccupancies[i] = np.sum(
#                RasterDict['RasterMatrix_b'][:,OccupiedBinsSet[i]], 
#                axis=0)
#        
#    (nTrials, nBins) = RasterDict['RasterMatrix_b'].shape
#    BinCountProportions = nBinOccupancies/nTrials
#    
#    return {
#            'OccupiedBinsSet' : OccupiedBinsSet,
#            'nBinOccupancies' : nBinOccupancies,
#            'BinCountProportions' : BinCountProportions,
#            'TimingEntropy' : - np.sum(BinCountProportions*np.log2(BinCountProportions))
#            }
    
def TimingEntropyCalculator(RasterDict):
    
    PostStimFilt = RasterDict['BinLeftEdgesVec'] >= 0.
    
    BinIndices = np.arange(0, RasterDict['BinLeftEdgesVec'][PostStimFilt].size)
    
    UniqueSequencesSet = np.unique(RasterDict['RasterMatrix_b'][:, BinIndices], axis=0)
    
    nSequenceInstances = np.nan*np.ones_like(UniqueSequencesSet[:,0])
    
    nSequenceIndices = np.arange(0, nSequenceInstances.size)
    
    for i in nSequenceIndices:
        
        DupArray = np.outer(np.ones_like(RasterDict['RasterMatrix_b'][:,0]), 
                            UniqueSequencesSet[i])
        
        nSequenceInstances[i] = np.sum(np.sum(
                RasterDict['RasterMatrix_b'] == DupArray, axis=1) == BinIndices.size, axis=0)
        
    SequenceProportions = nSequenceInstances/nSequenceInstances.sum()
    
    return {
            'Instances' : RasterDict['RasterMatrix_b'][:, BinIndices],
            'UniqueInstancesSet' : UniqueSequencesSet,
            'SummedInstances' : nSequenceInstances.sum(),
            'InstanceProportions' : SequenceProportions,
            'Entropy' : - np.sum(SequenceProportions*np.log2(SequenceProportions))
            }
    
def GenNestedDefaultDict():
    return defaultdict(GenNestedDefaultDict)


    
#def CountsMutInfDictAssembler(DataDict, BoundaryWindow, TargetBinWidth):
#    
#    # Cycle through the dictionary events and neurons fields.  Determine
#    # the probability of each stimulus (event type) being delivered.
#    
#    Resolution = 0.0005
#    
#    EventsList = list(DataDict['events'].keys())
#    eIndices = np.arange(0, len(EventsList))
#    
#    NeuronsList = list(DataDict['neurons'].keys())
#    nIndices = np.arange(0, len(NeuronsList))
#    
#    #OutDict = GenNestedDefaultDict()
#    CountsProportionsDicts = defaultdict(dict)
#    TimingProportionsDicts = defaultdict(dict)
#    
#    for n in nIndices:
#        
#        for e in eIndices:
#            
#            PeriEventTimesArray = PeriEventTimeStampExtractor(
#                    DataDict['events'][EventsList[e]], 
#                    DataDict['neurons'][NeuronsList[n]], 
#                    BoundaryWindow)
#   
#            RasterDict = RasterMatrixCompiler(PeriEventTimesArray, 
#                                              BoundaryWindow, 
#                                              Resolution)
#
#            RasterDict = RasterBinWidener(RasterDict, TargetBinWidth)
# 
#            CountEntropyDict = CountEntropyCalculator(RasterDict)
#
#
#            CountsProportionsDicts[NeuronsList[n]][EventsList[e]] = {
#                    'count_entropy' : CountEntropyDict['CountEntropy'],
#                    'proportions' : CountEntropyDict['CountProportions'],
#                    'nStimuli' : CountEntropyDict['SummedCounts']
#                    }
#    
#    return CountsProportionsDicts

def TimingMutInfDictAssembler(DataDict, BoundaryWindow, TargetBinWidth):

    # Cycle through the dictionary events and neurons fields.  Determine
    # the probability of each stimulus (event type) being delivered.
    
    Resolution = 0.0005
    
    EventsList = list(DataDict['events'].keys())
    eIndices = np.arange(0, len(EventsList))
    
    NeuronsList = list(DataDict['neurons'].keys())
    nIndices = np.arange(0, len(NeuronsList))
    
    #OutDict = GenNestedDefaultDict()
    TimingProportionsDicts = defaultdict(dict)
    
    CompositeSequenceSet = np.array([])
    
    for n in nIndices:
        
        for e in eIndices:
            
            PeriEventTimesArray = PeriEventTimeStampExtractor(
                    DataDict['events'][EventsList[e]], 
                    DataDict['neurons'][NeuronsList[n]], 
                    BoundaryWindow)
   
            RasterDict = RasterMatrixCompiler(PeriEventTimesArray, 
                                              BoundaryWindow, 
                                              Resolution)

            RasterDict = RasterBinWidener(RasterDict, TargetBinWidth)
 
            TimingEntropyDict = TimingEntropyCalculator(RasterDict)
            
            TimingProportionsDicts[NeuronsList[n]][EventsList[e]] = {
                    'timing_entropy' : TimingEntropyDict['Entropy'],
                    'UniqueInstancesSet' : TimingEntropyDict['UniqueInstancesSet'],
                    'proportions' : TimingEntropyDict['InstanceProportions'],
                    'nStimuli' : TimingEntropyDict['SummedInstances']
                    }
            
            if CompositeSequenceSet.size == 0:
                
                CompositeSequenceSet = TimingEntropyDict['UniqueInstancesSet']
                
            else: 
                
                CompositeSequenceSet = np.vstack((CompositeSequenceSet, 
                                              TimingEntropyDict['UniqueInstancesSet']))

    # Rerun through TimingProportionsDicts, concatenate all sequence sets into
    # a pool (above for-loop), keep only unique patterns 
    # (CompositeUniqueSequenceSet) compare the sequence set of each neuron
    # event pair to the composite set and find the indices in the composite 
    # set that point to each of the neuron-event sequence set.
        
    CompositeSequenceSet = np.unique(CompositeSequenceSet, axis=0)
    (nSequences, nBins) = CompositeSequenceSet.shape
    PointingIndices = np.arange(0, nSequences)
    OnesVec = np.ones_like(PointingIndices)
    
    for n in nIndices:
        
        for e in eIndices:
            
            ProportionsInGroupSet = np.zeros_like(PointingIndices).astype(float)
            
            r_s_pair_sequences = TimingProportionsDicts[
                    NeuronsList[n]][EventsList[e]]['UniqueInstancesSet']
            
            for i in np.arange(0, r_s_pair_sequences.shape[0]):
                
                DupArray = np.outer(OnesVec, r_s_pair_sequences[i])
                
                Filt = np.sum((CompositeSequenceSet == DupArray), 
                              axis=1) == r_s_pair_sequences.shape[1]
                
                ProportionsInGroupSet[PointingIndices[Filt]] = \
                        TimingProportionsDicts[NeuronsList[n]][
                                EventsList[e]]['proportions'][i].astype(float)
                        
            TimingProportionsDicts[
                    NeuronsList[n]][EventsList[e]]['group_proportions'] = \
                    ProportionsInGroupSet
        
    return TimingProportionsDicts

def CountsMutInfDictAssembler(DataDict, BoundaryWindow, TargetBinWidth):

    # Cycle through the dictionary events and neurons fields.  Determine
    # the probability of each stimulus (event type) being delivered.
    
    Resolution = 0.0005
    
    EventsList = list(DataDict['events'].keys())
    eIndices = np.arange(0, len(EventsList))
    
    NeuronsList = list(DataDict['neurons'].keys())
    nIndices = np.arange(0, len(NeuronsList))
    
    #OutDict = GenNestedDefaultDict()
    CountProportionsDicts = defaultdict(dict)
    
    CompositeSequenceSet = np.array([])
    
    for n in nIndices:
        
        for e in eIndices:
            
            PeriEventTimesArray = PeriEventTimeStampExtractor(
                    DataDict['events'][EventsList[e]], 
                    DataDict['neurons'][NeuronsList[n]], 
                    BoundaryWindow)
   
            RasterDict = RasterMatrixCompiler(PeriEventTimesArray, 
                                              BoundaryWindow, 
                                              Resolution)

            RasterDict = RasterBinWidener(RasterDict, TargetBinWidth)
 
            CountEntropyDict = CountEntropyCalculator(RasterDict)
            
            CountProportionsDicts[NeuronsList[n]][EventsList[e]] = {
                    'count_entropy' : CountEntropyDict['Entropy'],
                    'UniqueInstancesSet' : CountEntropyDict['UniqueInstancesSet'],
                    'proportions' : CountEntropyDict['InstanceProportions'],
                    'nStimuli' : CountEntropyDict['SummedInstances']
                    }
            
            if CompositeSequenceSet.size == 0:
                
                CompositeSequenceSet = CountEntropyDict['UniqueInstancesSet']
                
            else: 
                
                CompositeSequenceSet = np.vstack((CompositeSequenceSet, 
                                              CountEntropyDict['UniqueInstancesSet']))

    # Rerun through CountProportionsDicts, concatenate all sequence sets into
    # a pool (above for-loop), keep only unique patterns 
    # (CompositeUniqueSequenceSet) compare the sequence set of each neuron
    # event pair to the composite set and find the indices in the composite 
    # set that point to each of the neuron-event sequence set.
        
    CompositeSequenceSet = np.unique(CompositeSequenceSet, axis=0)
    (nSequences, nBins) = CompositeSequenceSet.shape
    PointingIndices = np.arange(0, nSequences)
    OnesVec = np.ones_like(PointingIndices)
    
    for n in nIndices:
        
        for e in eIndices:
            
            ProportionsInGroupSet = np.zeros_like(PointingIndices).astype(float)
            
            r_s_pair_sequences = CountProportionsDicts[
                    NeuronsList[n]][EventsList[e]]['UniqueInstancesSet']
            
            for i in np.arange(0, r_s_pair_sequences.shape[0]):
                
                DupArray = np.outer(OnesVec, r_s_pair_sequences[i])
                
                Filt = np.sum((CompositeSequenceSet == DupArray), 
                              axis=1) == r_s_pair_sequences.shape[1]
                
                ProportionsInGroupSet[PointingIndices[Filt]] = \
                        CountProportionsDicts[NeuronsList[n]][
                                EventsList[e]]['proportions'][i].astype(float)
                        
            CountProportionsDicts[
                    NeuronsList[n]][EventsList[e]]['group_proportions'] = \
                    ProportionsInGroupSet
        
    return CountProportionsDicts

def MutualInformationCalculator(r_given_s_Dict):
    
    StimuliTypesList = list(r_given_s_Dict.keys())
    nStimuliTypes = len(StimuliTypesList)
    
    StimuliTypesIndices = np.arange(0, nStimuliTypes)
    
    # Initialize P[s]
    P_of_s = np.nan*np.ones(nStimuliTypes)
    
    for i in StimuliTypesIndices:
    
        P_of_s[i] = r_given_s_Dict[StimuliTypesList[i]]['nStimuli']
        
    P_of_s = P_of_s / np.sum(P_of_s)
    
    # Calc P[r]
    
    P_of_r = np.zeros_like(r_given_s_Dict[StimuliTypesList[0]]['group_proportions'])

    for i in StimuliTypesIndices:
        
        P_of_r += P_of_s[i] * r_given_s_Dict[StimuliTypesList[i]]['group_proportions']
        
    I_m_array = np.zeros_like(P_of_r)
    
    for i in StimuliTypesIndices:
        
        P_of_r_given_s = r_given_s_Dict[StimuliTypesList[i]]['group_proportions']
        
        Filt = (P_of_r_given_s != 0.) 
        I_m_array[Filt] += P_of_s[i] * P_of_r_given_s[Filt] * np.log2(P_of_r_given_s[Filt] / P_of_r[Filt])
            
    return np.sum(I_m_array)

def JointMutInfDictAssembler(DataDict, BoundaryWindow, TargetBinWidth):
    # Cycle through the dictionary events and neurons fields.  Determine
    # the probability of each stimulus (event type) being delivered.
    
    Resolution = 0.0005    
    
    EventsList = list(DataDict['events'].keys())
    eIndices = np.arange(0, len(EventsList))
    
    NeuronsList = list(DataDict['neurons'].keys())
    nIndices = np.arange(0, len(NeuronsList))    
    
    JointCountProportionsDicts = defaultdict(dict)
    
    # Initialize the set to contain the group collection of instances from all events
    CompositeSequenceSet = np.array([])
    
    for e in eIndices:
        
        # Initialize sub dictionary to contain individual responses of each neuron
        CombinedInstances = np.array([])
        
        for n in nIndices:
            
            # Aquire instances about each event for each neuron.  Put in sub
            # directory to be combined outside this loop
            
            PeriEventTimesArray = PeriEventTimeStampExtractor(
                    DataDict['events'][EventsList[e]], 
                    DataDict['neurons'][NeuronsList[n]], 
                    BoundaryWindow)
   
            RasterDict = RasterMatrixCompiler(PeriEventTimesArray, 
                                              BoundaryWindow, 
                                              Resolution)

            RasterDict = RasterBinWidener(RasterDict, TargetBinWidth)
 
            CountEntropyDict = CountEntropyCalculator(RasterDict)
#            
#            CountInstancesDict[NeuronsList[n]] = {
#                    'Instances' : CountEntropyDict['Instances'],
#                    'count_entropy' : CountEntropyDict['Entropy'],
#                    'UniqueInstancesSet' : CountEntropyDict['UniqueInstancesSet'],
#                    'proportions' : CountEntropyDict['InstanceProportions'],
#                    'nStimuli' : CountEntropyDict['SummedInstances']
 #                   }
            
            if CombinedInstances.size == 0:
                
                CombinedInstances = CountEntropyDict['Instances']
                
            else:
                
                CombinedInstances = np.hstack((CombinedInstances, 
                                               CountEntropyDict['Instances']))
                
        # Make a mock raster dict that can be run through Timing Entropy
        # routine
        MockRasterDict = {
                'RasterMatrix_b' : CombinedInstances,
                'BinLeftEdgesVec' : np.arange(0, CombinedInstances.shape[1])
                }
        
        TimingEntropyDict = TimingEntropyCalculator(MockRasterDict)
        
        JointCountProportionsDicts[EventsList[e]] = {
                'timing_entropy' : TimingEntropyDict['Entropy'],
                'UniqueInstancesSet' : TimingEntropyDict['UniqueInstancesSet'],
                'proportions' : TimingEntropyDict['InstanceProportions'],
                'nStimuli' : TimingEntropyDict['SummedInstances']
                }
        
        if CompositeSequenceSet.size == 0:
            
            CompositeSequenceSet = TimingEntropyDict['UniqueInstancesSet']
            
        else: 
            
            CompositeSequenceSet = np.vstack((CompositeSequenceSet, 
                                          TimingEntropyDict['UniqueInstancesSet']))
            
    CompositeSequenceSet = np.unique(CompositeSequenceSet, axis=0)
    (nSequences, nBins) = CompositeSequenceSet.shape
    PointingIndices = np.arange(0, nSequences)
    OnesVec = np.ones_like(PointingIndices)    

    for e in eIndices:
        
        ProportionsInGroupSet = np.zeros_like(PointingIndices).astype(float)
        
        r_s_pair_sequences = JointCountProportionsDicts[EventsList[e]]['UniqueInstancesSet']
        
        for i in np.arange(0, r_s_pair_sequences.shape[0]):
            
            DupArray = np.outer(OnesVec, r_s_pair_sequences[i])
            
            Filt = np.sum((CompositeSequenceSet == DupArray), 
                          axis=1) == r_s_pair_sequences.shape[1]
            
            ProportionsInGroupSet[PointingIndices[Filt]] = \
                    JointCountProportionsDicts[
                            EventsList[e]]['proportions'][i].astype(float)
                    
        JointCountProportionsDicts[
                EventsList[e]]['group_proportions'] = \
                ProportionsInGroupSet
                               
#        JointInstancesByEventDicts[event] = MockRasterDict
                
    return JointCountProportionsDicts

############# Homework 3 functions ############
def zScoreInputFormatter(PeriEventActivityDict):
    
    KeysList = list(PeriEventActivityDict.keys())
    
    InputColumns =  np.empty((1, 1))

    for k in KeysList:

        InputColumn = np.array([PeriEventActivityDict[k].flatten()]).transpose()

        if InputColumns.shape == (1, 1):
            
            InputColumns = InputColumn
            
        else:
            
            InputColumns = np.hstack((InputColumns, InputColumn))
        
    return InputColumns

def zScoreInputFormatter2(PeriEventActivityDicts):
    
    OuterKeysList = list(PeriEventActivityDicts.keys())
        
    InnerKeysList = list(PeriEventActivityDicts[OuterKeysList[0]].keys())
    
    OuterInputColumns = np.empty((1, 1))
    
    for j in OuterKeysList:
    
        InnerInputColumns =  np.empty((1, 1))
        
        for k in InnerKeysList:
    
            InputColumn = np.array([PeriEventActivityDicts[j][k].flatten()]).transpose()
    
            if InnerInputColumns.shape == (1, 1):
                
                InnerInputColumns = InputColumn
                
            else:
                
                InnerInputColumns = np.hstack((InnerInputColumns, InputColumn))
        
        if OuterInputColumns.shape == (1, 1):
            
            OuterInputColumns = InnerInputColumns
            
        else:
            
            OuterInputColumns = np.vstack((OuterInputColumns, InnerInputColumns))
        
    return OuterInputColumns    
    
def zScoreColumns(InputColumns):
    
    MeansRow = np.mean(InputColumns, axis=0)
    StdDevRow = np.std(InputColumns, axis=0, ddof=1)
    OnesVec =  np.ones_like(InputColumns[:,0].transpose())
    
    zScoredInputColumns = \
        (InputColumns - np.outer(OnesVec, MeansRow)) / np.outer(OnesVec, StdDevRow) 
    
    return zScoredInputColumns

def OutputFormatter(OutputColumns, nBins):
    
    (OutputColumnLength, nNeurons) = OutputColumns.shape
    nTrials = int(OutputColumnLength/nBins)
    
    OutputArray = np.empty((nTrials, nBins*nNeurons))
    
    BinsSepVec = np.hstack([np.arange(0, OutputColumnLength,  nBins), OutputColumnLength])
    BinsSepVec_ind = np.arange(0, BinsSepVec.size)
    
    # Reformat data for covariance matrix (nTrials X nBins*nPCs)
    #for i in BinsSepVec_ind[0:-2]:
    for i in BinsSepVec_ind[0:-1]:
        
        RowVec = np.reshape(OutputColumns[BinsSepVec[i]:BinsSepVec[i + 1], :].transpose(), 
                            (nBins*nNeurons,), order='C')
        
        OutputArray[i, :] = RowVec
        
    return OutputArray

def ExtractPC_ProjectionPlots(PCA_PlotData, nTrialsPerEvent, nPCs, EventsList):
    
    # Cumulatively sum the nTrialsPerEvent array to get slicing boundaries of
    # the different event type sets.
    SliceBoundaries = np.int_(np.hstack((0, np.cumsum(nTrialsPerEvent))))
    
    # Generate slicing boundaries to extract trials corresponding each event type
    (nTotalTrials, nBinsTimesPCs) = PCA_PlotData.shape
    
    # Initialize the array to contain trial-averaged PC projection plots for
    # each event type
    PC_Plot_avgs = np.empty((SliceBoundaries.size - 1, nBinsTimesPCs))
    
    # Compile trial-averaged PC-projection plots
    for i in np.arange(0, (SliceBoundaries.size - 1)):
        
        PC_Plot_avgs[i,:] = np.sum(PCA_PlotData[
                SliceBoundaries[i]:SliceBoundaries[i + 1],:], 
                axis=0)/nTrialsPerEvent[i]
    
    # Generate slicing boundaries to extract each PC projection from the
    # concatenated average
    PC_width = nBinsTimesPCs/nPCs
    SliceBoundaries = np.int_(np.hstack((np.arange(0, nPCs*PC_width, PC_width), 
                                 nPCs*PC_width)))
    
    # Initialize the dictionary to contain each PC for each event type.
    PC_PlotsDict = defaultdict(dict)
    
    # Slice the PC projection averages and write to dictionary.
    for i in np.arange(0, len(EventsList)):
    
        for j in np.arange(0,(SliceBoundaries.size - 1)):
            
            PC_PlotsDict[EventsList[i]]['pc_'+str(j + 1)+'_psth'] = \
                PC_Plot_avgs[i, SliceBoundaries[j]:SliceBoundaries[j + 1]]
                
    return PC_PlotsDict

############# Homework 4 functions ############
    
def EuclideanDistanceCalculator(vec1, vec2):
    
    #return np.linalg.norm(vec2 - vec1)
    DiffArray = (vec2 - vec1)
    return np.sqrt(np.sum(DiffArray**2))

def MakeTemplate(PEA_Array):
    
    (nTrials, nBinsByNeurons) = PEA_Array.shape
    
    return np.sum(PEA_Array, axis=0)/nTrials

def LeaveOneOutParser(PEA_Array, EventSetIndices, TestTrialIndex):
    
    Filt = (EventSetIndices != TestTrialIndex)
    
    TemplateAverage = MakeTemplate(PEA_Array[EventSetIndices[Filt], :])
    
    return PEA_Array[TestTrialIndex, :], TemplateAverage

def MutualInformationFromConfusionMatrix(ConfusionMatrix):
    
    # Compute mutual information
    P_Joint = ConfusionMatrix/np.sum(np.sum(ConfusionMatrix, axis=1), axis=0)
    
    P_True = np.sum(P_Joint, axis=1)
    
    P_Predicted = np.sum(P_Joint, axis=0)
    
    MutInf = 0.
    
    for i in np.arange(0, P_True.size):
    
        for j in np.arange(0, P_Predicted.size):
            
            if (P_True[i] != 0.) & (P_Predicted[j] != 0.) & (P_Joint[i, j] != 0.):
                
                MutInf += P_Joint[i, j]*np.log2(P_Joint[i, j] / 
                             (P_True[i] * P_Predicted[j]))
                
    return MutInf

def PerformanceCalculator(ConfusionMatrix):
    
    nTotalTrials = np.sum(np.sum(ConfusionMatrix, axis=1), axis=0)
    
    nCorrectTrials = np.trace(ConfusionMatrix)
    
    return nCorrectTrials / nTotalTrials

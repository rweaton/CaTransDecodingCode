#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 12:56:48 2019

@author: thugwithyoyo
"""

from scipy.stats import t, spearmanr
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

import PeriEventTraceFuncLib as PETFL

#nBins = 5
#x = DifferenceMagnitudes[0,:]
#y = DifferenceMagnitudes[1,:]

def MutInfCalculator(x, y, nBins):

    # Not sure if this approach is correct.  Need to confirm before using.    
    (JointProbHist, xedges, yedges) = np.histogram2d(x, y, bins=nBins)
    
    JointProbHist = JointProbHist / np.sum(np.sum(JointProbHist))
    
    (NumXBins, NumYBins) = JointProbHist.shape
    
    MargHist_X = np.sum(JointProbHist, axis=1)
    
    MargHist_Y = np.sum(JointProbHist, axis=0)
    
    MutInf = 0.
    
    xIndices = np.arange(0, NumXBins)
    yIndices = np.arange(0, NumYBins)
    
    xFilt = (MargHist_X != 0.)
    yFilt = (MargHist_Y != 0.)
    
    for x in xIndices[xFilt]:
        
        for y in yIndices[yFilt]:
            
            if (JointProbHist[x, y] > 0.):
                
                Scale = (xedges[x+1] - xedges[x])*(yedges[y+1] - yedges[y])
                
                MutInf += Scale*JointProbHist[x, y]*np.log2(JointProbHist[x, y]/
                                       (MargHist_X[x]*MargHist_Y[y]))
            
            else:
                
                pass
            
    return MutInf

#MutInf = MutInfCalculator(x, y, nBins)
    
def TuningAndMappingParser(CellFluorTraces_Frame, Tuning_Frame, 
                           CentroidFrame, ParamsDict):
    
    ##### Get constants
    # Set whether or not to eliminate non-tuned cells (eliminate cells with tuning
    # in middle quantiles)

    
    # Acquire cell names from column headers
    CellIDs = CellFluorTraces_Frame.columns.values[1:]
    
    # Count total number of cells in recording.
    (NumCells,) = CellIDs.shape
      
    RemovalFilt = np.full((NumCells,), False, dtype=bool)
    
    Filt = np.logical_or(np.isnan(CentroidFrame[CellIDs].loc['x_coord'].values), 
                         np.isnan(CentroidFrame[CellIDs].loc['y_coord'].values))
    
    RemovalFilt = np.logical_or(RemovalFilt, Filt)
    

    if ParamsDict['RemoveNontunedCells'] == True:
    
        # Filter by tuning values; remove cells with tuning within middle quantile range
#        InMiddleQuantiles = ((Tuning_Frame['ScalarTuningIndex'].values > ParamsDict['DistQuantiles'][0]) & 
#                             (Tuning_Frame['ScalarTuningIndex'].values < ParamsDict['DistQuantiles'][1]))
#        FilteredCellIDs = CellIDs[~InMiddleQuantiles]
        
        for ntd in ParamsDict['NontunedDomains']:
            
            Filt = ((Tuning_Frame['ScalarTuningIndex'].values > ntd[0]) &
                    (Tuning_Frame['ScalarTuningIndex'].values < ntd[1]))
            
            RemovalFilt = np.logical_or(RemovalFilt, Filt)
    
    FilteredCellIDs = CellIDs[~RemovalFilt]
    
    # Generate a unique list of pairs from entries in CellIDs
    UniqueIndexPairsList, UniquePairsList = PETFL.UniquePairsGenerator(FilteredCellIDs)
    
    # Count number of pairs
    (NumPairs, _) = UniquePairsList.shape
    
    # Initialize array to contain magnitudes of tuning differences and
    # differences in cell location vectors.
    DifferenceMagnitudes = np.empty((2, NumPairs))
     
    # Compute magnitudes of tuning differences for the cell pairs
    DifferenceMagnitudes[1, :] = np.abs(
            Tuning_Frame['ScalarTuningIndex'].loc[UniquePairsList[:,1]].values - 
            Tuning_Frame['ScalarTuningIndex'].loc[UniquePairsList[:,0]].values)
    
    # Compute magnitudes of vector differences in position between cell pairs
    DifferenceMagnitudes[0, :] = np.linalg.norm(
            CentroidFrame[UniquePairsList[:,1]].loc[['x_coord', 'y_coord']].values - 
            CentroidFrame[UniquePairsList[:,0]].loc[['x_coord', 'y_coord']].values,
            axis=0)
    
#    MaxVals = np.array([np.max(DifferenceMagnitudes, axis=1)]).transpose()
    
#    MaxValsArray = np.dot(MaxVals, np.ones((1,NumPairs)))
    
#    NormalizedDifferenceMagnitudes = np.divide(DifferenceMagnitudes, MaxValsArray)
    
    MapIndicesFromTuningSort = \
        np.argsort(Tuning_Frame['ScalarTuningIndex'].loc[FilteredCellIDs].values)
    
    SortedFilteredCellIDs = FilteredCellIDs[MapIndicesFromTuningSort]

    SortedTuningVals = Tuning_Frame['ScalarTuningIndex'].loc[
            SortedFilteredCellIDs].values
            
    SortedCentroids = CentroidFrame[SortedFilteredCellIDs].loc[
            ['x_coord', 'y_coord']].values
            
    (NumFilteredCells,) = SortedFilteredCellIDs.shape
    
    IndexSet = np.arange(0, NumFilteredCells)

    IndexSetsByTuningDict = defaultdict()
    
    for domain in ParamsDict['TuningDomains']:
        
        IndexSetsByTuningDict[str(domain)] = IndexSet[
                (SortedTuningVals >= domain[0]) &
                (SortedTuningVals < domain[1])]
    
    return {
            'DifferenceMagnitudes': DifferenceMagnitudes,
            'SortedFilteredCellIDs': SortedFilteredCellIDs,
            'SortedTuningVals': SortedTuningVals,
            'SortedCentroids': SortedCentroids,
            'IndexSetsByTuningDict': IndexSetsByTuningDict
            }
    
def CalcDistSeparatingGroupCenters(SortedCentroids, IndexSetsByTuningDict):
    
    (NumDims, NumCentroids) = SortedCentroids.shape
    NumSets = len(IndexSetsByTuningDict)
    Centers = np.empty((NumDims, NumSets))
    
    i = 0
    for IndexSet in IndexSetsByTuningDict:
        
        (NumIndices,) = IndexSetsByTuningDict[IndexSet].shape
        Centers[:,i] = np.sum(SortedCentroids[:,IndexSetsByTuningDict[IndexSet]], 
                               axis=-1)/NumIndices
        i += 1
        
    Centers = np.transpose(Centers)
    
    CentersMat = np.empty((NumSets, NumSets, NumDims))
    OnesVec = np.ones((1, NumSets))
    
    for j in np.arange(0, NumDims):
    
        CentersMat[:, :, j] = Centers[:, j:(j+1)] @  OnesVec
    
    DiffMat = CentersMat - np.transpose(CentersMat, axes=(1,0,2))
    
    return {
            'GroupCenters': Centers,
            'GroupSeparationDistances': np.linalg.norm(DiffMat, axis=2)
           }
    
def ComputeLinearCorrelationStats(x, y):
    
    # Make sure that the sample data are in the form of numpy arrays
    x = np.array(x)
    y = np.array(y)
    
    # Initialize output dictionary.
    StatsDict = defaultdict()
    
    # Make sure that input argument arrays have same number of elements
   
    # Count number of (paired) entries; the number of datapoints defined by 
    (StatsDict['NumPairs'],) = x.shape
    
    try:
        
        # Compute covariance matrix
        StatsDict['CovMat'] = np.cov(x, y)
        
        # Compute correlation matrix
        StatsDict['CorrCoefMat'] = np.corrcoef(x, y)
    
    except:
        
        print('Error: Unable to generate covariance or correlation matrix.\nLikely input arrays are different sizes.')

    # Calculate linear regression fit coefficients
    beta1 = StatsDict['CovMat'][0,1] / StatsDict['CovMat'][0,0]
    beta0 = np.mean(y) - beta1*np.mean(x)
    StatsDict['LinFitCoefs'] = np.array([beta0, beta1])
    
    # Extract Pearson product-moment correlation coefficient.
    StatsDict['PearsonCoef'] = StatsDict['CorrCoefMat'][0,1]
    
    # Calculate the corresponding t-statistic for the correlation coefficient 
    StatsDict['Pearson_t'] = StatsDict['PearsonCoef'] * (
            np.sqrt(StatsDict['NumPairs'] - 2) / 
            np.sqrt(1. - StatsDict['PearsonCoef']**2))
    
    # Calculate the degrees of freedom
    StatsDict['DegOfFreedom'] = StatsDict['NumPairs'] - 2
    
    # Retrieve p-value from two-tailed t-distribution
    StatsDict['Pearson_pval'] = t.sf(np.abs(StatsDict['Pearson_t']), 
             StatsDict['DegOfFreedom'])*2  # two-sided pvalue = Prob(abs(t)>tt)
    
    return StatsDict

def ComputeNonLinearCorrelationStats(x, y):
    
    # Make sure that the sample data are in the form of numpy arrays
    x = np.array(x)
    y = np.array(y)
    
    # Initialize output dictionary.
    StatsDict = defaultdict()
    
    # Make sure that input argument arrays have same number of elements
   
    # Count number of (paired) entries; the number of datapoints defined by 
    (StatsDict['NumPairs'],) = x.shape

    # Calculate Spearman rank-order correlation coefficient
    [StatsDict['SpearmanCoef'], StatsDict['Spearman_pval']] = spearmanr(x, y)
    
    # Calculate the degrees of freedom
    StatsDict['DegOfFreedom'] = StatsDict['NumPairs'] - 2
    
    return StatsDict
    
def CorrelationShuffledBootstrap(x, y, NumReps, CorrType):
        
    # Make sure that the sample data are in the form of numpy arrays
    x = np.array(x)
    y = np.array(y)
    
    # Initialize array to hold bootstrap simulation output
    Corr_hat = np.empty((NumReps,), dtype=float)
    
    if CorrType == 'Pearson':

        for r in np.arange(0, NumReps):  
            
            # Shuffle order of values in y array to destroy correspondence
            # between entries of x and y arrays
            np.random.shuffle(y)
            
            # Calculate and retain the Pearson product-moment correlation 
            # coefficient using the two arrays (x and y) with broken 
            # correspondence.
            Corr_hat[r] = ComputeLinearCorrelationStats(x, y)['PearsonCoef']
            
    if CorrType == 'Spearman':
        
        for r in np.arange(0, NumReps):
            
            # Shuffle order of values in y array to destroy correspondence
            # between entries of x and y arrays
            np.random.shuffle(y)
            
            # Calculate and retain the Pearson product-moment correlation 
            # coefficient using the two arrays (x and y) with broken 
            # correspondence.
            Corr_hat[r] = ComputeNonLinearCorrelationStats(x, y)['SpearmanCoef']            
            
        #Corr_hat = np.hstack((Corr_hat, NonLinStatsDict['SpearmanCoef']))
        
    return Corr_hat

def SeparationShuffledBootstrap(SortedCentroids, IndexSetsByTuningDict, NumReps):
    
    (NumDims, NumCentroids) = SortedCentroids.shape
    NumSets = len(IndexSetsByTuningDict)
    
    Sep_hat = np.empty((NumReps, NumSets, NumSets))
    
    CentroidsSelectionIndices = np.arange(0, NumCentroids)
    
    for i in np.arange(0, NumReps):
        
        np.random.shuffle(CentroidsSelectionIndices)
        
        GroupSeparationsDict = CalcDistSeparatingGroupCenters(
                SortedCentroids[:, CentroidsSelectionIndices], IndexSetsByTuningDict)
        
        Sep_hat[i, :, :] = GroupSeparationsDict['GroupSeparationDistances']
        
    return Sep_hat
        
def ChanceNullHypothTest(ObservedVal, SampleDist, Alpha, TailType):
    
    #StatsDict = defaultdict()

    SampleDist = np.hstack((SampleDist, ObservedVal))
    
    NumReps = SampleDist.shape[0]
    
    if TailType == 'double':
        
        QuantileLevels = np.array([Alpha/2., 1. - Alpha/2.])
        CountFilt = (np.abs(SampleDist) >= ObservedVal)
        
    elif TailType == 'low':
        
        QuantileLevels = np.array([Alpha])
        CountFilt = (SampleDist < ObservedVal)
        
    elif TailType == 'high':
        
        QuantileLevels = np.array([1. - Alpha])
        CountFilt = (SampleDist > ObservedVal)
        
    Bootstrap_pval = (np.sum(CountFilt)/NumReps)
    RejectNullHypoth = (Bootstrap_pval < Alpha)
    
    return {
            'NumReps': NumReps,
            'Bootstrap_pval': Bootstrap_pval,
            'QuantileLevels': QuantileLevels,
            'RejectNullHypoth': RejectNullHypoth
            }

def ScatterPlotGenerator(x, y, PlotParams, **kwargs):
    
    # Initialize figure and axes object.  If one is provided as one of the
    # optional keyword arguments, use it.  If not, generate new figure and
    # axes object.
    if kwargs['AxesHandle']:
        
        axs = kwargs['AxesHandle']
        plt.sca(axs)
        fig = plt.gcf()
        
    else:
        
        fig, axs = plt.subplots(1,1)

    # Generate scatter plot
    axs.scatter(x, y,
                marker=PlotParams['MarkerType'],
                s=PlotParams['MarkerSize'], 
                alpha=PlotParams['MarkerOpacity'], 
                c=PlotParams['MarkerColor'],
                label=PlotParams['Labels'][PlotParams['GroupIndex']])
    
    # Remove bounding box
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    # Generate labels for x and y axes
    axs.set_ylabel(PlotParams['yLabel'])
    axs.set_xlabel(PlotParams['xLabel'])
    axs.set_aspect(PlotParams['AspectRatio'])
    
    TextBoxDict = defaultdict() 
    
    if 'LinCorrStats' in kwargs:
        
        x_new = np.linspace(np.min(x), np.max(x), 50)
        y_new = (x_new*kwargs['LinCorrStats']['LinFitCoefs'][1] + 
                 kwargs['LinCorrStats']['LinFitCoefs'][0])
    
        axs.plot(x_new, y_new, 'k--')

        TextBoxDict['NumPairs'] = kwargs['LinCorrStats']['NumPairs']    
        TextBoxDict['PearsonCoef'] = kwargs['LinCorrStats']['PearsonCoef']
        TextBoxDict['Pearson_pval'] = kwargs['LinCorrStats']['Pearson_pval']
        

            
    if 'NonLinCorrStats' in kwargs:
        TextBoxDict['NumPairs'] = kwargs['NonLinCorrStats']['NumPairs']
        TextBoxDict['SpearmanCoef'] = kwargs['NonLinCorrStats']['SpearmanCoef']
        TextBoxDict['Spearman_pval'] = kwargs['NonLinCorrStats']['Spearman_pval']

        
    if 'GroupSeparationsDict' in kwargs:
        
        CenterVec = kwargs['GroupSeparationsDict']['GroupCenters'][PlotParams['GroupIndex']]
        axs.scatter(CenterVec[0], CenterVec[1],
#                    marker=PlotParams['MarkerType'],
                    marker='+',
                    s=4*PlotParams['MarkerSize'], 
                    alpha=PlotParams['MarkerOpacity'], 
                    c=PlotParams['MarkerColor'],
                    label=PlotParams['Labels'][PlotParams['GroupIndex']] + '_centroid')
        
 #       LegendText['SepDistance'] = kwargs['GroupSeparationsDict'][
 #               'GroupSeparationDistances'][]]
    
    # Assemble info for axes text box
    TextBoxString = ''
    
    for d in TextBoxDict:
        
        if d[-5:] == '_pval':
            TextBoxString = TextBoxString + d + '={:.2E}'.format(TextBoxDict[d]) + '\n'
            
        elif d[-4:] == 'Coef':
            TextBoxString = TextBoxString + d + '={:0.4f}'.format(TextBoxDict[d]) + '\n'
            
        elif d[:3] == 'Num':
            TextBoxString = TextBoxString + d + '={:5d}'.format(TextBoxDict[d]) + '\n'
        
        else:
            TextBoxString = TextBoxString + d + '={:6.5f}'.format(TextBoxDict[d]) + '\n'
    
    # build a rectangle in axes coords
    left, width = .25, .5
    bottom, height = .25, .5
    
    # Report statistics from correlation analysis.
    right = left + width
    top = bottom + height
    axs.text(left, top, TextBoxString, transform=axs.transAxes, fontsize=8.0)
    
    return fig, axs
    
def SampleDistribPlotter(ObservedVal, SampleDist, NumBins, Alpha, **kwargs):
    
    # Append observed value to sample distribution
    SampleDist = np.hstack((SampleDist, ObservedVal))
    
    # Recount the number of samples comprising the sample distribution
    (NumSamples,) = SampleDist.shape
    
    # Initialize figure and axes object.  If one is provided as one of the
    # optional keyword arguments, use it.  If not, generate new figure and
    # axes object.
    if kwargs['AxesHandle']:
        
        axs = kwargs['AxesHandle']
        plt.sca(axs)
        fig = plt.gcf()
        
    else:
        
        fig, axs = plt.subplots(1,1)
    
    # Generate histogram of sampling distribution.  Included title and labels
    # for x and y axes.
    Counts, BinEdges = np.histogram(SampleDist, bins=NumBins)
    BarWidth = BinEdges[1] - BinEdges[0]
    axs.bar(BinEdges[0:-1], Counts/NumSamples, align='edge', width=BarWidth)
    axs.set_title('Bootstrapped Sample Distribution')
    #axs.set_ylabel('probability density')
    axs.set_xlabel('sample value')
    
    # Remove bounding box
    axs.spines['right'].set_visible(False)
    axs.spines['top'].set_visible(False)
    
    # Get limits of y-axis of current plot
    yMax = axs.get_ylim()[1]
    
    # Compute the median of the sample distribution and plot as vertical 
    #dashed line
    SampleDist_median = np.median(SampleDist)
    axs.plot(np.array([SampleDist_median, SampleDist_median]), 
             np.array([0., yMax]),'k--', label='median')
        
    # Plot the observed value of the statistic as a dash-dotted red vertical
    # line on over the sample distribution.
    axs.plot(np.array([ObservedVal, ObservedVal]), np.array([0., yMax]),
             'r-.', label='observed')
    
    # build a rectangle in axes coords
    left, width = .25, .5
    bottom, height = .25, .5
    right = left + width
    top = bottom + height

    if 'HypothTestDict' in kwargs:
        
        # Retrieve quantile levels and compute quantile boundaries 
        #QuantileLevels = np.array([Alpha/2., 1. - Alpha/2.])
        Quantiles = np.quantile(SampleDist, kwargs['HypothTestDict']['QuantileLevels'])
        
        # Plot as vertical dotted lines, each of the quanitles calculated above.
        for q in np.arange(0, Quantiles.shape[0]):
            
            Label = str(kwargs['HypothTestDict']['QuantileLevels'][q]) + ' quant.'
            axs.plot(np.array([Quantiles[q], Quantiles[q]]), np.array([0., yMax]),
                     'k:', label=Label)
        
        # Report values from hypotheses test in a text box.
        axs.text(right, top, 'Obs.=%4.3f\nR=%6d\np-val.=%.2E' % 
                 (ObservedVal,
                  kwargs['HypothTestDict']['NumReps'],
                  kwargs['HypothTestDict']['Bootstrap_pval']), 
                 transform=axs.transAxes,
                 fontsize=8.0
                 )
        
    else:
        
        # Report only the number of samples (repetitions) in the bootstrapped
        # sampling distribution.
        axs.text(right, top, 'Obs.=%4.3f\nR=%6d' % (ObservedVal, NumReps), 
                 transform=axs.transAxes,
                 fontsize=8.0)

    axs.legend(prop={'size': 6})
    
    return fig, axs

#################################
#### Begin script execution #####
#################################    
#TuningDomainBound = np.min(np.abs(DistQuantiles))

ParsingParamsDict = {
         'RemoveNontunedCells': False,
         'DistQuantiles': DistQuantiles,
#         'TuningDomains': np.array([[-np.inf, DistQuantiles[0]],
#                                    [DistQuantiles[0], DistQuantiles[1]],
#                                    [DistQuantiles[1], np.inf]])
#         'TuningDomains': np.array([[-np.inf, -TuningDomainBound],
#                                    [-TuningDomainBound, TuningDomainBound],
#                                    [TuningDomainBound, np.inf]])
         'TuningDomains': np.array([[-np.inf, -ParamsDict['TuningCutoffLevel']],
                                    [-ParamsDict['TuningCutoffLevel'], ParamsDict['TuningCutoffLevel']],
                                    [ParamsDict['TuningCutoffLevel'], np.inf]]),
         'NontunedDomains': np.array([[-ParamsDict['TuningCutoffLevel'], 0.],
                                     [0., ParamsDict['TuningCutoffLevel']]])
         }
            
ParserDict = TuningAndMappingParser(CellFluorTraces_Frame, Tuning_Frame, 
                           CentroidFrame, ParsingParamsDict)

TuningGroupNames = list(ParserDict['IndexSetsByTuningDict'].keys())

# Plot cell position in FOV scatter plot and bootstrap analysis on
# distance separating group center of mass positioning

GroupSeparationsDict = CalcDistSeparatingGroupCenters(
                            ParserDict['SortedCentroids'], 
                            ParserDict['IndexSetsByTuningDict'])

GroupSepsMat = GroupSeparationsDict['GroupSeparationDistances']

#GroupsToAnalyze = np.array([TuningGroupNames[0], TuningGroupNames[2]])
#(NumGroupsToAnalyze,) = GroupsToAnalyze.shape
#for g in GroupsToAnalyze:


NumReps = 4999

Sep_hat = SeparationShuffledBootstrap(ParserDict['SortedCentroids'], 
                                      ParserDict['IndexSetsByTuningDict'], 
                                      NumReps)

Alpha = 0.05
TailType = 'high'
SepHypothTestDict = ChanceNullHypothTest(GroupSepsMat[0, 2], Sep_hat[:, 0, 2], 
                                         Alpha, TailType)

NumBins = 50
fig, axs = plt.subplots(1,2)
SampleDistribPlotter(GroupSepsMat[0, 2], Sep_hat[:, 0, 2], 
                  NumBins, Alpha, HypothTestDict=SepHypothTestDict,
                  AxesHandle=axs[1])

PlotParams = {
              'yLabel':'vertical coord. (pixels)',
              'xLabel': 'horizontal coord. (pixels)',
              'AspectRatio': 1.,
              'MarkerType': 'o',
              'MarkerSize': 10.,
              'MarkerOpacity': 0.5,
              'MarkerColor': 'lightgray',
              'Labels': TuningGroupNames,
              'GroupIndex': 1
             }

CellLocs = ParserDict['SortedCentroids'][:,
                     ParserDict['IndexSetsByTuningDict'][
                             TuningGroupNames[PlotParams['GroupIndex']]]]

fig, axs[0] = ScatterPlotGenerator(CellLocs[0,:], CellLocs[1,:], PlotParams, 
                AxesHandle=axs[0], GroupSeparationsDict=GroupSeparationsDict)

PlotParams['MarkerColor'] = 'magenta'
PlotParams['GroupIndex'] = 2

CellLocs = ParserDict['SortedCentroids'][:,
                     ParserDict['IndexSetsByTuningDict'][
                             TuningGroupNames[PlotParams['GroupIndex']]]]

fig, axs[0] = ScatterPlotGenerator(CellLocs[0,:], CellLocs[1,:], PlotParams, 
                AxesHandle=axs[0], GroupSeparationsDict=GroupSeparationsDict)

PlotParams['MarkerColor'] = 'green'
PlotParams['GroupIndex'] = 0

CellLocs = ParserDict['SortedCentroids'][:,
                     ParserDict['IndexSetsByTuningDict'][
                             TuningGroupNames[PlotParams['GroupIndex']]]]

fig, axs[0] = ScatterPlotGenerator(CellLocs[0,:], CellLocs[1,:], PlotParams, 
                AxesHandle=axs[0], GroupSeparationsDict=GroupSeparationsDict)

axs[1].set_xlabel('distance separating group centers (pixels)')
fig.suptitle(File[0:19] + ' Cell FOV position and target preference')
axs[0].set_title('Cell position in FOV')

if ParsingParamsDict['RemoveNontunedCells']:

    FilteredString = '_Filtered'

else:

    FilteredString = ''
    
fig.savefig(Path+os.sep + File[0:19] + '_TuningGroupCoMSeparation_' + 
             ParamsDict['TuningType']+ FilteredString + '.svg')

# End of cell location scatter plots

# Begin correlation analysis (and plot) between distance separating cell pairs
# in FOV and difference in target tuning of cell pairs. 
LinStatsDict = ComputeLinearCorrelationStats(ParserDict['DifferenceMagnitudes'][0,:], 
                                             ParserDict['DifferenceMagnitudes'][1,:])

NonLinStatsDict = ComputeNonLinearCorrelationStats(ParserDict['DifferenceMagnitudes'][0,:], 
                                                   ParserDict['DifferenceMagnitudes'][1,:])

NumReps = 4999
Corr_hat = CorrelationShuffledBootstrap(ParserDict['DifferenceMagnitudes'][0,:], 
                                        ParserDict['DifferenceMagnitudes'][1,:], 
                                        NumReps, 'Spearman')

Alpha = 0.05
TailType = 'double'
HypothTestDict = ChanceNullHypothTest(NonLinStatsDict['SpearmanCoef'], 
                                      Corr_hat, Alpha, TailType)

NumBins = 50
Alpha = 0.05
fig, axs = plt.subplots(1,2)
fig, axs[1] = SampleDistribPlotter(NonLinStatsDict['SpearmanCoef'], Corr_hat, 
                            NumBins, Alpha, HypothTestDict=HypothTestDict,
                            AxesHandle=axs[1])

PlotParams = {
              'xLabel': r'|$\vec r_{cell 2} - \vec r_{cell 1}$| (pixels)',
              'yLabel': r'|$TI_{cell 2} - TI_{cell 1}$|',
              'AspectRatio': 'auto',
              'MarkerType': 'o',
              'MarkerSize': 1.,
              'MarkerOpacity': 0.5,
              'MarkerColor': 'gray',
              'Labels': np.array(['differences between cell pair']),
              'GroupIndex': 0
             }

fig, axs[0] = ScatterPlotGenerator(ParserDict['DifferenceMagnitudes'][0,:], 
                                   ParserDict['DifferenceMagnitudes'][1,:], 
                                   PlotParams, 
                                   AxesHandle=axs[0], 
                                   LinCorrStats=LinStatsDict, 
                                   NonLinCorrStats=NonLinStatsDict)

axs[1].set_xlabel('Spearman corr. coef.')
axs[1].set_title('Bootstrapped sample distribution', fontsize=10)
axs[0].set_title('Tuning vs. spatial\ndifferences among cell pairs', fontsize=10)
fig.suptitle(File[0:19] + ' Spatial/Tuning Correlation')

if ParsingParamsDict['RemoveNontunedCells']:
    
    FilteredString = '_Filtered'

else:
    
    FilteredString = ''
    
fig.savefig(Path+os.sep + File[0:19] + '_Tuning-SpatialCorrelation_' + 
             ParamsDict['TuningType']+ FilteredString + '.svg')
# End of correlation scatter plot analysis
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 07:15:43 2019

@author: thugwithyoyo
"""

import PeriEventTraceFuncLib
from scipy.stats import t

# Set constants
FitNormalized = False

# Set whether or not to eliminate non-tuned cells (eliminate cells with tuning
# in middle quantiles)
RemoveUntunedCells = True

CellIDs = CellFluorTraces_Frame.columns.values[1:]

# Count total number of cells in recording.
(NumCells,) = CellIDs.shape

# Filter by tuning values; remove cells with tuning within middle quantile range
InMiddleQuantiles = ((Tuning_Frame['ScalarTuningIndex'].values > DistQuantiles[0]) & 
                     (Tuning_Frame['ScalarTuningIndex'].values < DistQuantiles[1]))

if RemoveUntunedCells:

    FilteredCellIDs = CellIDs[~InMiddleQuantiles]
    
else:
    
    FilteredCellIDs = CellIDs
    
# Generate a unique list of pairs from entries in CellIDs
UniqueIndexPairsList, UniquePairsList = UniquePairsGenerator(FilteredCellIDs)

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

MaxVals = np.array([np.max(DifferenceMagnitudes, axis=1)]).transpose()

MaxValsArray = np.dot(MaxVals, np.ones((1,NumPairs)))

NormalizedDifferenceMagnitudes = np.divide(DifferenceMagnitudes, MaxValsArray)

fig, axs = plt.subplots(nrows=1,ncols=1)

# Compute correlation between the two measures. User selects above whether to 
# plot and fit the normalized or the raw (non-normalized) magnitudes.
if FitNormalized == True:

    axs.scatter(NormalizedDifferenceMagnitudes[0,:], 
                NormalizedDifferenceMagnitudes[1,:],s=0.5)
    
    (z, residuals, rank, singular_values, rcond) = np.polyfit(
                   NormalizedDifferenceMagnitudes[0,:], 
                   NormalizedDifferenceMagnitudes[1,:], 1, 
                   full=True)
    
    x_new = np.linspace(np.min(NormalizedDifferenceMagnitudes[0,:]),
                        np.max(NormalizedDifferenceMagnitudes[0,:]), 
                        50)
    
    CorrCoefMat = np.corrcoef(NormalizedDifferenceMagnitudes[0,:],
                              NormalizedDifferenceMagnitudes[1,:])
    
    # Compute the t-test p-value for the Pearson correlation coefficient
    # Retrieve Pearson product-moment correlation coefficient
    r = CorrCoefMat[0,1]
    
    # Calculate the corresponding t-statistic for the correlation coefficient 
    tt = r * np.sqrt(NumPairs - 2) / np.sqrt(1. - r**2)
    
    # Calculate the degrees of freedom
    df = NumPairs - 2
    
    # Retrieve p-value from two-tailed t-distribution
    pval = t.sf(np.abs(tt), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print( 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval))
    
else:

    axs.scatter(DifferenceMagnitudes[0,:], 
                DifferenceMagnitudes[1,:],s=0.5)
    
    (z, residuals, rank, singular_values, rcond) = np.polyfit(
                   DifferenceMagnitudes[0,:], 
                   DifferenceMagnitudes[1,:], 1, 
                   full=True)
    
    x_new = np.linspace(np.min(DifferenceMagnitudes[0,:]),
                        np.max(DifferenceMagnitudes[0,:]), 
                        50)
    
    CorrCoefMat = np.corrcoef(DifferenceMagnitudes[0,:],
                              DifferenceMagnitudes[1,:])   
    
    # Compute the t-test p-value for the Pearson correlation coefficient
    # Retrieve Pearson product-moment correlation coefficient
    r = CorrCoefMat[0,1]
    
    # Calculate the corresponding t-statistic for the correlation coefficient 
    tt = r * np.sqrt(NumPairs - 2) / np.sqrt(1. - r**2)
    
    # Calculate the degrees of freedom
    df = NumPairs - 2
    
    # Retrieve p-value from two-tailed t-distribution
    pval = t.sf(np.abs(tt), df)*2  # two-sided pvalue = Prob(abs(t)>tt)
    print( 't-statistic = %6.3f pvalue = %6.4f' % (tt, pval))
    
f = np.poly1d(z)

axs.plot(x_new, f(x_new), 'k--')

# build a rectangle in axes coords
left, width = .25, .5
bottom, height = .25, .5

# Report statistics from correlation analysis.
right = left + width
top = bottom + height
axs.text(right, top, 'corr.=%6.3f\np-val.=%6.5f\n$n_{pairs}$=%6d' % 
         (r, pval, NumPairs), transform=axs.transAxes)
axs.set_xlabel(r'|$\vec r_{cell 2} - \vec r_{cell 1}$| (pixels)')
axs.set_ylabel(r'|$TI_{cell 2} - TI_{cell 1}$| (a.u.)')
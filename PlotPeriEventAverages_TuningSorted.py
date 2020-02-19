#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 23:40:34 2019

@author: thugwithyoyo
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from collections import defaultdict


import PeriEventTraceFuncLib as PETFL
#from CaTraceNormalizer import *
import CaTraceNormalizer as CTN

import os
import tkinter as tk
from tkinter.filedialog import askopenfilename

#RestoreFilePath = SavePath +'.dat'
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()

# Determine parent directory and filename from complete path.
Drive, Path_and_file = os.path.splitdrive(RestoreFilePath)
Path, File = os.path.split(Path_and_file)

#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-18-11-20-21/SW_400ms_100ms/2018-12-18-11-20-21_new_unique_400ms_SamFiltered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-18-11-20-21_new_unique_400ms_SamFiltered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-17-11-38-42_new_unique_400ms_SamFiltered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-19-11-24-58_new_unique_400ms_SamFiltered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56_unique_400ms_SamFiltered.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56/SW_400ms_100ms/2018-12-10-11-37-56_unique_400ms_SlidingWindow.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-14-11-01-41/JointHemisphereDecoding/SingleArmDirectionTask/2018-12-14-11-01-41_new_unique_SW_LHem.dat'
#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-14-11-01-41_new_unique_400ms_SamFiltered.dat'
try:
    
    exec(open('./RestoreShelvedWorkspaceScript.py').read())

except:
    
    print('Error when restoring shelved workspace.  Will attempt to continue')
    
### Specify a dictionary of parameters for tuning calcluation and plotting.
# Time domain to consider for calculation of tuning index
ParamsDict['SearchDomain'] = [-1., 1.]

# Specify the method to be used to calculate tuning method. 
#ParamsDict['TuningType'] = 'PlainDifference'
#ParamsDict['TuningType'] = 'WeightedDifference'
ParamsDict['TuningType'] = 'DiffOfAvgsOverSumOfMagOfAvgs'

# For the WeightedDifference method, specify the standard deviation of the 
# Gaussian weighting function.
ParamsDict['StdDev'] = 0.75

# Specify the number of bins to use to generate the histogram.
ParamsDict['NumHistoBins'] = 20

# Specify the list of quantiles to plot as vertical lines on histogram
ParamsDict['QuantilesList'] = [0.25, 0.75]

# The magnitude of the tuning index value that separates tuned cells from 
# non-tuned cells.
ParamsDict['TuningCutoffLevel'] = 0.75
#ParamsDict['TuningCutoffLevel'] = 0.5
#ParamsDict['TuningCutoffLevel'] = 0.2

# Pie chart color coding
ParamsDict['PieColors'] = np.array(['green', 'lightgray', 'magenta'])

# Name of the color scheme used for plotting heat map.
ColorMap = 'seismic'

# 2018-12-11-10-53-04 cells to select
#CellsToSelectList = ['C086', 'C103', 'C017', 'C053', 'C074']

# 2019-01-07-10-31-41 cells to select
#CellsToSelectList = ['C07', 'C29', 'C43', 'C69', 'C28', 'C36', 'C59', 
#                     'C04', 'C21', 'C54', 'C59']

# 2019-01-07-10-45-52 cells to select
#CellsToSelectList = ['C64', 'C01', 'C29', 'C39', 'C27', 'C28', 'C30']

# 2018-12-18-11-20-21 cells to select
#CellsToSelectList = ['C031', 'C085', 'C025', 'C081', 'C099', 'C101']
#CellsToSelectList = ['C085', 'C025', 'C031']

# 2018-12-20-11-26-21 cells to select
#CellsToSelectList = ['C085', 'C023', 'C073']

#2018-12-20-11-26-34 cells to select
#CellsToSelectList = ['C07', 'C03', 'C09']

RefEventsList = RefEventsDict['RefEventsList']
SortingRefEvent = RefEventsList[1]

nEvents = len(RefEventsList)

# Load in CentroidFrame just for the 2018-12-11-... analysis. Shelve is wonkey.
# !!!! Be sure to comment this out after!!!
#CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)

# Acquire number of columns in fluorescence trace dataframe
nCols = CellFluorTraces_Frame.shape[1]

# zScore all calcium traces prior to peri-event extraction
FullTraces = pd.DataFrame(CellFluorTraces_Frame._slice(slice(1, nCols), 1)).values.transpose()

zScoredTraces = CTN.zScoreTraces(FullTraces, ParamsDict).transpose()

# Replace raw fluorescence traces with z-scores in the trace dataframe
IndicesOfColsToChange = np.arange(1, nCols)
ColNames = list(CellFluorTraces_Frame.columns)

for i in IndicesOfColsToChange:
    CellFluorTraces_Frame[ColNames[i]] = zScoredTraces[:, i - 1]
    
#CellFluorTraces_Frame._slice(slice(1, nCols), 1).values = zScoredTraces.transpose() 

# Recompile PeriEvent activity dict to include Full Domain in boundary window.
PeriEventExtractorDict = PETFL.PeriEventExtractor_Trace(BehavDict, 
                                                  CellFluorTraces_Frame, 
                                                  RefEventsDict, 
                                                  ParamsDict['BoundaryWindow'])

# Count the number of trials and determine length of each row in the formatted
# Peri-event activity array.
NumTrials, NumCellsByNumSamples = PeriEventExtractorDict['PEA_Array'].shape

# Count the total number of columns in the fluorescence trace dataframe
(_, NumColumns) = CellFluorTraces_Frame.shape

# Cell traces are columns after the first timestamp column at index 0
NumCells = NumColumns - 1

# Determine the number of entries within each peri-event snippet
NumSamples = int(NumCellsByNumSamples/NumCells)

# Generate the x-axis event relative time vector for plotting
RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
                         ParamsDict['BoundaryWindow'][1], 
                         num=NumSamples, endpoint=False)

# Count the number of event types to be plotted
NumEvents = len(RefEventsList)

# Initialize the three dimensional array to contain traces to  average
AveragedTracesMatrices = np.empty((NumCells, NumSamples, NumEvents))
#PeakNormalizedTraces = np.empty((NumCells, NumSamples, NumEvents))

# Initialize dicts to contain information for sorting traces
SortProcessingDict = defaultdict(dict)
AbsSortProcessingDict = defaultdict(dict)
NormVals = np.zeros((NumCells,))

###############################################################################
################### BEGIN SORTING AVERAGED Z-SCORED TRACES ####################
###############################################################################
# Iterate through each event type in list
# Set index counter to initial value
i = 0
for RefEvent in RefEventsList:
    
    # Generate filter for selecting trials of the current event type
    TracesFilt = PeriEventExtractorDict['TrialIndicesByEventDict'][RefEvent]
    
    # Average together trials of the current event type
    AveragedTraces = np.mean(PeriEventExtractorDict['PEA_Array'][TracesFilt,:],
                                      axis=0)
    
    # Reshape flattend average into 2D array and write to the collection be
    # plotted.
    AveragedTracesMatrices[:,:,i] = np.reshape(AveragedTraces, (NumCells, NumSamples))
    
    
#    SortProcessingDict[RefEvent] = CaTraceSorter(AveragedTracesMatrices[:,:,i],
#                                                 ParamsDict, 'PeakSort')
    
#    PeakNormalizedTraces[:, :, i] = PeakNormalizer(AveragedTracesMatrices[:,:,i])

#    SortProcessingDict[RefEvent] = CaTraceSorter(PeakNormalizedTraces[:,:,i],
#                                                 ParamsDict, 'ActivitySort')

    # Sort trace by absolute value.  (NOT SURE IF THIS STILL NEEDS TO BE PERFORMED???)
    AbsSortProcessingDict[RefEvent] = CTN.CaTraceSorter(
            AveragedTracesMatrices[:,:,i], ParamsDict, 'AbsPeakSort')
    
    NormVals = NormVals + np.abs(AbsSortProcessingDict[RefEvent]['MaxVals'])
    
    # Increment index to select next event type on next iteration
    i += 1

NormVals = NormVals / len(RefEventsList)
    
#vmin = np.min(np.min(np.min(AveragedTracesMatrices, axis=0),axis=0),axis=0)

# Determine maximum value to set color map.  Take floor to the next 0.5 level
vmax = np.max(np.max(np.max(AveragedTracesMatrices, axis=0),axis=0),axis=0)
vmax = np.ceil(2.0*vmax)/2.0
vmin = -vmax

#vmin = np.min(np.min(np.min(PeakNormalizedTraces, axis=0),axis=0),axis=0)
#vmax = np.max(np.max(np.max(PeakNormalizedTraces, axis=0),axis=0),axis=0)

#TuningMatrix = np.divide(np.diff(AveragedTracesMatrices, axis=-1)[:,:,0],
#                         np.dot(np.array([NormVals]).transpose(), 
#                         np.array([np.ones(NumSamples,)])))

# Begin calculation of scalar tuning indices for z-score averaged traces from 
# each cell
#(NumTraces, NumSamples) = AveragedTracesMatrices[:,:,0].shape
#
#RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
#                         ParamsDict['BoundaryWindow'][1], 
#                         num=NumSamples, endpoint=False)
#
#SearchDomainFilt = (RelTimeVec >= ParamsDict['SearchDomain'][0]) & \
#                   (RelTimeVec <= ParamsDict['SearchDomain'][1])

# Calculate scalar tuning values for all averaged traces according to the method 
# specified in ParamsDict['TuningType']
TuningScalars = CTN.GetTuning(AveragedTracesMatrices[:,:,::-1], 
                              ParamsDict['TuningType'], 
                              ParamsDict=ParamsDict)
     
#TuningScalars = GetTuning(AveragedTracesMatrices[:,:,::-1], 'DiffOfAvgOverSumOfAvgMag', ParamsDict=ParamsDict)
#TuningMatrix = GetTuning(AveragedTracesMatrices[:,:,::-1], 'MaxAbsDiffMatrixNorm')[:,:,0]
#SortProcessingDict = CaTraceSorter(TuningMatrix, ParamsDict, 'PeakSort')
#SortProcessingDict = CaTraceSorter(TuningMatrix, ParamsDict, 'AreaSort')
#SortProcessingDict = CaTraceSorter(TuningMatrix, ParamsDict, 'AvgSort')

# Generate the index map to use to to traces by tuning index 
SortProcessingDict['AvgVals'] = np.array([TuningScalars]).transpose()
SortProcessingDict['SortIndices'] = np.argsort(SortProcessingDict['AvgVals'], axis=0).transpose()[0]

# Write cell identities, tuning values and traces ordering in dataframe.

#CellLabels = np.array(list(CellFluorTraces_Frame.columns.values))[1:]
#Tuning_Frame = pd.DataFrame(index=CellLabels, columns=['ScalarTuningIndex', 'HeatMapRowIndex'])
#Tuning_Frame['ScalarTuningIndex'] = SortProcessingDict['AvgVals'].transpose()[0]
#
## Perform reverse mapping operation to indicate each the row index of the 
## associated trace in the tuning-sorted heat map.
#IndexList = np.arange(0,NumCells)
#for i in IndexList:
#    
#    Filt = (SortProcessingDict['SortIndices'] == i) 
#    
#    Tuning_Frame.at[CellLabels[i], 'HeatMapRowIndex'] = IndexList[Filt][0]
    
# Above commented out section put in a subroutine in CaTraceNormalizer.
# Generate a dataframe of tuning values for each cell included in the 
# fluorescence dataframe. 
Tuning_Frame = CTN.TuningByCellFrameGen(CellFluorTraces_Frame, 
                                        SortProcessingDict, 
                                        ParamsDict)

# Save tuning dataframe into Excel file.
Tuning_Frame.to_excel(Path+os.sep+File[0:19]+'_TIsByCell_' + 
                      str(ParamsDict['SearchDomain']) + '_' +
                      ParamsDict['TuningType']+'.xlsx')

###############################################################################
########################## BEGIN COLORMAPS PLOTTING ###########################
###############################################################################
# Initialize figure and constituent subplot axes
fig, axes  = plt.subplots(1, NumEvents)

# generate two 2D grids within x & y bounds for plotting heat map data
xv, yv = np.meshgrid(RelTimeVec, np.arange(1, NumCells + 1, 1))

# Iterate over all event types in list and generate a heat map for each in its
# own subplot

# Set starting index to reference array of subplot pointers
i = 0
for RefEvent in RefEventsList:
    
    im = axes[i].pcolor(xv, yv, AveragedTracesMatrices[SortProcessingDict['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)    
#    im = axes[i].pcolor(xv, yv, AveragedTracesMatrices[SortProcessingDict[SortingRefEvent]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)
#    im = axes[i].pcolormesh(xv, yv, PeakNormalizedTraces[SortProcessingDict[RefEventsList[0]]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap='coolwarm')

    # select cells
    CellsToSelectList = np.array(CellsToSelectList)
    (NumCellsToSelect,) = CellsToSelectList.shape
    
    for j in np.arange(0, NumCellsToSelect):
        
        # Set the y-axis value to plot the dotted lines.  NOTE: Friggin matplotlib
        # starts the colormap plot at index 1, not 0, so each of these have to 
        # be incremented by 1 to have the correct correspondence.
        ylevel_low = Tuning_Frame['HeatMapRowIndex'][CellsToSelectList[j]] + 1
        ylevel_high = ylevel_low + 2
        
        axes[i].plot([RelTimeVec[0], RelTimeVec[-1]],[ylevel_low, ylevel_low], 
            linewidth=0.5, linestyle='dashed', label=CellsToSelectList[j],
            color='purple')
        
        #axes[i].plot([RelTimeVec[0], RelTimeVec[-1]],[ylevel_high, ylevel_high], 
        #    linewidth=0.5, linestyle='dashed', label=CellsToSelectList[j])
        
        axes[i].text(RelTimeVec[0]-0.5, ylevel_low, CellsToSelectList[j], 
            fontsize=6, verticalalignment='center')
        
    #fig.colorbar(im, ax=axes[i])
    axes[i].set_xlabel('time relative to reach (sec.)')
    
    if i == 0:
        axes[i].set_ylabel('order index')
        
    if i > 0:
        
        axes[i].set_yticks([])

#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M7T1_Entry_ts']
#ParamsDict['RefEventsList'] = ['M6T0_Entry_ts', 'M6T1_Entry_ts']    
    if (RefEvent == 'M6T0_Entry_ts'):
        SubPlotTitle = 'RH_Zone1'
    elif(RefEvent == 'M6T1_Entry_ts'):
        SubPlotTitle = 'RH_Zone2'
    elif(RefEvent == 'M7T0_Entry_ts'):
        SubPlotTitle = 'LH_Zone1'
    elif(RefEvent == 'M7T1_Entry_ts'):
        SubPlotTitle = 'LH_Zone2'
        
    axes[i].title.set_text(SubPlotTitle)
    
    i += 1

# Assign colors across relevant numerical  range.
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
#num = int(((vmax - vmin)/0.5) + 1.)

vmin = np.floor(vmin)/2.0
bounds = np.linspace(vmin,  vmax, num=301, endpoint=True)

# Plot color bar key for heat maps
cb = fig.colorbar(im, ax=list(axes), orientation='vertical', boundaries=bounds, 
                  cmap=ColorMap, norm=norm)
# Still need to get rid of bounding box, change grid labels to increment to 
# values of 0.5, and thin out the rectangle surrounding the bar...

#cb = mpl.colorbar.ColorbarBase(
#        fig.axes[2], boundaries=np.linspace(vmin,  vmax, num=301, endpoint=True),
#        orientation='vertical', cmap=ColorMap, norm=norm, )

#
#mpl.colorbar.ColorbarBase(
#        cb, boundaries=np.linspace(vmin,  vmax, num=300, endpoint=False), 
#        orientation='horizontal', cmap=ColorMap, norm=norm)
Ticks = np.arange(vmin, vmax + 0.5, 0.5)
cb.set_ticks(Ticks)
cb.set_ticklabels(Ticks)
cb.set_label('average z-score')

# Define figure name
FigureTitle = File[0:19] + ' trial-averaged z-scores of $\Delta$F $Ca^{2+}$ response\n' + 'tuning time domain = ' + str(ParamsDict['SearchDomain']) + '\n'
fig.suptitle(FigureTitle)

# Set size of entire figure
fig.set_size_inches(11., 8.5)

# Save figure
fig.savefig(Path+os.sep+File[0:19] + '_TraceHeatMapsByZone_' + 
            str(ParamsDict['SearchDomain']) + '_' +
            ParamsDict['TuningType'] + '.svg')

#  Generate difference plot. 
#fig2, axs2 = plt.subplots(1,1)
#
#DiffMat = np.diff(AveragedTracesMatrices, axis=-1)[:,:,0]
##DiffMat = np.diff(PeakNormalizedTraces, axis=-1)[:,:,0]
#
#im2 = axs2.pcolormesh(xv, yv, DiffMat[SortProcessingDict[SortingRefEvent]['SortIndices'],:], cmap=ColorMap)
#fig2.colorbar(im2, ax=axs2)
#axs2.set_ylabel('Cell identity')
#axs2.set_xlabel('time relative to reach (sec.)')


##vmax = np.max(np.max(TuningMatrix, axis=0),axis=0)
#vmax = GetMax(TuningMatrix)
#vmax = np.floor(2.0*vmax)/2.0
#vmin = -vmax
#
## Generate a tuning plot
#fig2, axs2 = plt.subplots(1,1)
#im2 = axs2.pcolormesh(xv, yv, TuningMatrix[SortProcessingDict['SortIndices'],:], cmap='PiYG', vmin=vmin, vmax=vmax)
#cb2 = fig2.colorbar(im2, ax=axs2)
#cb2.set_label('tuning index')
#axs2.set_ylabel('Cell identity')
#axs2.set_xlabel('time relative to reach (sec.)')
#
#FigureTitle = file[0:19] + ' reach tuning by cell (z-scores of $Ca^{2+}$ $\Delta$F)'
#fig2.suptitle(FigureTitle)

###############################################################################
############################# BEGIN HISTOGRAM PLOTTING ########################
###############################################################################
# Define figure name.  This is defined below.
#FigureTitle = File[0:19] + ' distribution of cell tuning indices \nfrom peri-reach $Ca^{2+}$ response traces'
#FigureTitle = File[0:19] + ' proportions of cell tuning indices \nfrom peri-reach $Ca^{2+}$ response traces'

# Generate histogram of scalar tuning values.
# A figure with two subplot axes will illustrate the histogram on the left and
# the pie chart on the right.
fig2, axs2 = plt.subplots(nrows=1, ncols=2)
fig2.suptitle(FigureTitle)
BinBounds = np.linspace(-1.,1., num=(ParamsDict['NumHistoBins'] + 1), endpoint=True)
FrequencyCount, bins = np.histogram(SortProcessingDict['AvgVals'].transpose()[0], 
                                    bins=BinBounds)

# Compute positions along x axis of median and quantiles lines.
DistMedian = np.median(SortProcessingDict['AvgVals'])
DistQuantiles = np.quantile(SortProcessingDict['AvgVals'], ParamsDict['QuantilesList'])

# Generate bar graph of proportions of cell tuning 
Proportions = FrequencyCount/np.sum(FrequencyCount)
axs2[0].bar(BinBounds[0:-1], Proportions, align='edge', 
         width=(BinBounds[1] - BinBounds[0]),
         color='grey', edgecolor='dimgrey')

# Set y-axis boundary for histogram plot.
#ymax = int(np.ceil(1.1*np.max(FrequencyCount)))
ymax = np.ceil(10.*np.max(Proportions))/10.
ymin = 0

# Plot median line
axs2[0].plot([DistMedian, DistMedian], [0., ymax], color='black', 
          linestyle='dashed', linewidth=2.)

# Plot quantile(s) line(s)
(NumQuantiles,) = DistQuantiles.shape

for i in np.arange(0,NumQuantiles):
    
    axs2[0].plot([DistQuantiles[i], DistQuantiles[i]], [0., ymax],  
              color='black', linestyle='dotted', linewidth=2.)

# Set figure properties.
axs2[0].set_ylim(ymin=ymin, ymax=ymax)
#axs2[0].set_yticks(list(range(ymin, ymax + 1, 2)))
axs2[0].spines['top'].set_visible(False)
axs2[0].spines['right'].set_visible(False)
#axs2[0].set_xticks(bins)
axs2[0].set_xlabel('tuning index')
axs2[0].set_ylabel('proportion of T.I.s')

# Save figure
#fig2.savefig(Path+os.sep + File[0:19] + '_TuningIndexHisto_' + 
#             str(ParamsDict['NumHistoBins']) + 'Bins_' + 
#             ParamsDict['TuningCutoffLevel'] + '_TIBound_' +
#             ParamsDict['TuningType'] + '.svg')

###############################################################################
########################### BEGIN PIE-CHART PLOTTING ##########################
###############################################################################
# Set axis on which to plot the pie chart.  Pie chart to be plotted on the 
# right subplot axis of the figure initialized above.
#fig3, axs3 = plt.subplots(1,1)
axs3 = axs2[1]

# Initialize array to contain proportions for the three tuning categories: 
# (zone 1-tuned, zone 2-tuned and untuned).
Proportions = np.empty((3,))
(TotalNumIndices,_) = SortProcessingDict['AvgVals'].shape

# Calulate the proportion of Zone2-tuned cells
Proportions[0] = np.sum(SortProcessingDict['AvgVals'] <=
                        -1.*ParamsDict['TuningCutoffLevel'])/TotalNumIndices

# Calculate the proportion of Zone1-tuned cells
Proportions[2] = np.sum(SortProcessingDict['AvgVals'] >=  
                        ParamsDict['TuningCutoffLevel'])/TotalNumIndices

# Calculate the remaining proportion  of non-selective cells.
Proportions[1] = 1. - (Proportions[0] + Proportions[2])

# Generate the corresponding label for each entry in the list of proportions
Labels = ['TI $\leq$ '+str(-1.*ParamsDict['TuningCutoffLevel']),
#          str(-1.*ParamsDict['TuningCutoffLevel'])+' < TI < '+ str(ParamsDict['TuningCutoffLevel']),
          'TI < |' + str(ParamsDict['TuningCutoffLevel'])+'|',
          'TI $\geq$ '+str(ParamsDict['TuningCutoffLevel'])]

# Generate pie chart
axs3.pie(Proportions, labels=Labels, autopct='%1.1f%%', colors=ParamsDict['PieColors'])

# Post total number of cells, N, for the proportions
axs3.text(1., 0., 'N = '+str(TotalNumIndices), fontsize=12,  transform=axs3.transAxes)

#fig3.savefig(Path+os.sep + File[0:19] + '_TuningIndexProportionsPie_' + 
#             ParamsDict['TuningCutoffLevel'] + '_TIBound_' +
#             ParamsDict['TuningType']+'.svg')

###############################################################################
############################ FINISH UP FIGURE AND SAVE ########################
###############################################################################
# Post title of entire figure
FigureTitle = (File[0:19] + ' distribution and proportions of trace-based cell tuning indices\n' 
    + ' tuning time domain = ' + str(ParamsDict['SearchDomain']) 
#    + ', threshold for tuned cells = |%2.2f|' % (ParamsDict['TuningCutoffLevel'],)
    )

# Post title of entire figure
fig2.suptitle(FigureTitle)

# Set figure size
fig2.set_size_inches(11., 8.5)

# Save figure
fig2.savefig(Path+os.sep + File[0:19] + '_TuningIndexHistoAndPie_' +
             str(ParamsDict['SearchDomain']) + '_' +
             str(ParamsDict['TuningCutoffLevel']) + '_' +
             ParamsDict['TuningType'] + '.svg')

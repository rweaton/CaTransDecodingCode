#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 21:15:13 2020

@author: thugwithyoyo
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

###########
# Here: run PlotPeriEventAverages_TuningSorted.py to load left hem. workspace.
###########

CellsToSelectList_LH = CellsToSelectList

Tuning_Frame_LH = Tuning_Frame

SortProcessingDict_LH = SortProcessingDict

AveragedTracesMatrices_LH = AveragedTracesMatrices
                      
###########
# Here: run PlotPeriEventAverages_TuningSorted.py to load right hem. workspace.
###########

CellsToSelectList_RH = CellsToSelectList

Tuning_Frame_RH = Tuning_Frame

SortProcessingDict_RH = SortProcessingDict

AveragedTracesMatrices_RH = AveragedTracesMatrices

NumCells_RH = NumCells 

# Begin joining operations
RowIndexOffset = np.max(Tuning_Frame_LH['HeatMapRowIndex'].values) + 1

OffsetIndices = Tuning_Frame_RH['HeatMapRowIndex'].values + RowIndexOffset

Tuning_Frame_RH['HeatMapRowIndex'] = OffsetIndices

Tuning_Frame = pd.concat([Tuning_Frame_LH, Tuning_Frame_RH])
CellsToSelectList = np.hstack([CellsToSelectList_LH, CellsToSelectList_RH])

JoinedTracesMatrices = np.vstack([AveragedTracesMatrices_LH[SortProcessingDict_LH['SortIndices'],:,:], 
                                  AveragedTracesMatrices_RH[SortProcessingDict_RH['SortIndices'],:,:]])

NumCells = NumCells_LH + NumCells_RH

# Determine maximum value to set color map.  Take floor to the next 0.5 level
vmax = np.max(np.max(np.max(JoinedTracesMatrices, axis=0),axis=0),axis=0)
vmax = np.ceil(2.0*vmax)/2.0
vmin = -vmax


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
    
    im = axes[i].pcolor(xv, yv, JoinedTracesMatrices[:,:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)
#    im = axes[i].pcolor(xv, yv, AveragedTracesMatrices[SortProcessingDict['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)    
#    im = axes[i].pcolor(xv, yv, AveragedTracesMatrices[SortProcessingDict[SortingRefEvent]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap=ColorMap)
#    im = axes[i].pcolormesh(xv, yv, PeakNormalizedTraces[SortProcessingDict[RefEventsList[0]]['SortIndices'],:,i], vmin=vmin, vmax=vmax, cmap='coolwarm')

    # select cells
    CellsToSelectList = np.array(CellsToSelectList)
    (NumCellsToSelect,) = CellsToSelectList.shape
    
    axes[i].plot([RelTimeVec[0], RelTimeVec[0]], [RowIndexOffset, NumCells], 
            linewidth=8.0, linestyle='solid', label='Right hemisphere',
            color='cyan')
        
    axes[i].plot([RelTimeVec[0], RelTimeVec[0]], [1, RowIndexOffset], 
            linewidth=8.0, linestyle='solid', label='Left hemisphere',
            color='lightgray')
    
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

axes[0].set_ylabel('')
axes[0].set_yticks([])

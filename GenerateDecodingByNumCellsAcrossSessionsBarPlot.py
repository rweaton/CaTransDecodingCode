#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 18:35:05 2019

@author: thugwithyoyo
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import tkinter as tk
from tkinter.filedialog import askopenfilename

# This script generates a bar plot and line plot (with errorbars) of 
# decoding performance statistics from selected sessions. Workspaces selected
# must be those generated from the script: 
# DecodingPerformanceAcrossNumOfCells.py.

#### Set constants #####
# Integer width to position clusters of bars in bar plot.
ClusterIncrement = 1.

# The proportion of space to use to separate bar groups.
GroupSepWidthProportion = 0.2

# Set default colormap to determine color sequence as bars are plotted. Not
# sure if this is actually working...
mpl.rc('image', cmap='nipy_spectral')

# Specify the "number of included cells" datapoints to extract from performance
# arrays from each session.  Specifically, the particular filter below extracts
# "no. of cells included" values at 10, 20, 30, 40, 50 and all cells from the
# particular workspaces being opened.  This definately can and should be
# generalized for any selected workspace that they user happens to select, via
# name comparison and index filtering....
EntriesFilt = np.array([0, 5, 7, 9, 11, 13, -1])
EntriesFilt = EntriesFilt[::-1]

#### Initialize variables ####
# Initialize arrays to contain aggregate decoding performance statistics.
AgrPerfMeans = np.array([])
AgrPerfSEs = np.array([])
AgrShuffledPerfMeans = np.array([])
AgrShuffledPerfSEs= np.array([])

# Initialize array to contain list of names of selected sessions to plot.
SessionNames = np.array([])


# Initialize the termination case for subsequent session selection loop.
GetAnotherSession = True

# Ask user to navigate and select workspaces to be included.  Use tkinter
# gui dialog for user to select files and from which path info can be
# determined.
while (GetAnotherSession == True):
    
    # Acquire path of workspace to load.
    root = tk.Tk()
    RestoreFilePath = askopenfilename()
    root.withdraw()
    
    # Open workspace using Shelve loading script.  Place in try/except block
    # to try and bypass occasional loading errors.
    #exec(open('./RestoreShelvedWorkspaceScript.py').read())
    try:
    
        exec(open('./RestoreShelvedWorkspaceScript.py').read())
        
    except:
    
        print('Unshelving error.  Will attempt to continue...')
    
    # Determine parent directory and filename from complete path.
    drive, path_and_file = os.path.splitdrive(RestoreFilePath)
    path, file = os.path.split(path_and_file)
    
    # Grow composite arrays with each iteration of the loop.
    if AgrPerfMeans.shape[0] == 0:
        
        AgrPerfMeans = PerfMeans[EntriesFilt]
        AgrPerfSEs  = PerfSEs[EntriesFilt]
        AgrShuffledPerfMeans = ShuffledPerfMeans[EntriesFilt]
        AgrShuffledPerfSEs = ShuffledPerfSEs[EntriesFilt]
        
        SessionNames = file[0:19]
        
    elif AgrPerfMeans.shape[0] > 0:
        
        AgrPerfMeans = np.vstack([AgrPerfMeans, PerfMeans[EntriesFilt]])
        AgrPerfSEs  = np.vstack([AgrPerfSEs, PerfSEs[EntriesFilt]])
        AgrShuffledPerfMeans = np.vstack([AgrShuffledPerfMeans, ShuffledPerfMeans[EntriesFilt]])
        AgrShuffledPerfSEs = np.vstack([AgrShuffledPerfSEs, ShuffledPerfSEs[EntriesFilt]])
    
        SessionNames = np.hstack([SessionNames, file[0:19]])
        
    # Query user if another session is to be added to the session pool to be 
    # plotted.
    GetAnotherSession = tk.messagebox.askyesno(message='Include another session?')
    
    root.withdraw()

#### Begin bar plot figure generation ####
# Append to arrays performance statistics from the "shuffled performance" group.
BarNames = np.hstack([np.array(np.array(X[EntriesFilt], dtype=int), dtype=str), 'Shuffled'])
BarVals = np.hstack([AgrPerfMeans, np.array([AgrShuffledPerfMeans[:, 0]]).transpose()])
SEVals = np.hstack([AgrPerfSEs, np.array([AgrShuffledPerfSEs[:, 0]]).transpose()])

ShuffledBarVals = np.hstack([AgrShuffledPerfMeans, np.array([AgrShuffledPerfMeans[:, 0]]).transpose()])
ShuffledSEVals = np.hstack([AgrShuffledPerfSEs, np.array([AgrShuffledPerfSEs[:, 0]]).transpose()])

# Rename the first element of the bar sequence to be plotted.
BarNames[0] = 'All'

# Determine session and bars per group counts for subsequent iterative plot
# generation.
(NumSessions, NumBarsPerGroup) = BarVals.shape

# Calculate bar widths
BarWidth = ClusterIncrement*(1 - GroupSepWidthProportion)/NumBarsPerGroup

# Generate group position vector along x axis.
GroupPositions = np.arange(ClusterIncrement, (NumSessions + 1)*ClusterIncrement, 
                           ClusterIncrement)

BarIndices = np.arange(0, NumBarsPerGroup)

RelBarPositions = BarWidth*BarIndices - ClusterIncrement*(1 - GroupSepWidthProportion)/2.
RelErrorBarPositions = RelBarPositions + (BarWidth/2.)*np.ones_like(RelBarPositions)

fig1, ax = plt.subplots(nrows=1, ncols=1)

colormap = plt.get_cmap('nipy_spectral')
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=['black', 'blue', 'purple', 'red', 'orange', 'green', 'cyan']) 

for i in BarIndices:
    
    ax.bar(GroupPositions + RelBarPositions[i], BarVals[:,i], label=BarNames[i], align='edge', width=BarWidth)
    ax.errorbar(GroupPositions + RelErrorBarPositions[i], BarVals[:,i], yerr=SEVals[:,i], fmt='none', ecolor='black')

# Set plot axes limits, ticknames, legend location, etc...
ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
ax.set_xlabel('Session')
ax.set_ylabel('Decoding Accuracy')
ax.set_ylim([0., 1.])
NiceSessionNames = SessionNames
for i in np.arange(0, SessionNames.shape[0]):
    NiceSessionNames[i] = SessionNames[i][0:10]
    
ax.set_xticks(range(1, GroupPositions.shape[0]+1, 1))
ax.set_xticklabels(NiceSessionNames, rotation=45, ha='center')

# Remove plot "boundary box".
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Display axes title.
ax.set_title('Average trace decoding accuracy by no. of cells included across sessions')

#### Begin line plot of session averages across number of included cells plot. ####
MarkerSize='12'

fig2, ax2 = plt.subplots(nrows=1, ncols=1)

# Determine number of entries that comprise each session pool.
(NumEntries, _) = BarVals[:,1:-1].shape



# Calculate Standard Errors of the Means statistics to plot extents of errorbars
ySEs = np.divide(np.std(BarVals[:,1:-1], axis=0), np.sqrt(NumEntries))

#ax2.errorbar(np.arange(1, NumBarsPerGroup-1), np.mean(BarVals[:,1:-1], 
#             axis=0), yerr=np.std(BarVals[:,1:-1], axis=0), 
#             color='blue', ecolor='blue', fmt='.-', 
#             label='Decoding with subsampled units')
ErrorBarFormat = ''
# Plot errobars for the subsampled sets.
ax2.errorbar(np.arange(1, NumBarsPerGroup-1), np.mean(BarVals[:,1:-1], axis=0), 
             yerr=ySEs, color='orange', ecolor='orange', fmt=ErrorBarFormat,
             markersize=MarkerSize)

# Plot errorbar for the all cells included set.
NumEntries = BarVals[:,0].shape
ySEs = np.divide(np.std(BarVals[:,0], axis=0), np.sqrt(NumEntries))


#ax2.errorbar(np.array([0]), np.mean(BarVals[:,0], axis=0), 
#             yerr=np.std(BarVals[:,0], axis=0), color='black', 
#             ecolor='black', fmt='.-', label='Decoding with all units')
ax2.errorbar(np.array([0]), np.mean(BarVals[:,0], axis=0), 
             yerr=ySEs, color='orange', ecolor='orange', fmt=ErrorBarFormat,
             markersize=MarkerSize, label='Observed outcomes')

# Determine the value along the y-axis at which to plot the mean of 
#shuffled outcomes performance from sessions. 

#ShuffledBarVals = np.mean(AgrShuffledPerfMeans[:, EntriesFilt], axis=0)

# Calculate Standard Errors of the Shuffled Means statistics for extents of errorbars
#ShuffledySEs = np.divide(np.std(AgrShuffledPerfMeans[:, EntriesFilt], axis=0, ddof=1), np.sqrt(NumEntries))

#ax2.plot(np.array([0, NumBarsPerGroup-2]), np.array([ShuffledVal, ShuffledVal]), 'c--', label='Shuffled performance (all units)')


ySEs = np.divide(np.std(BarVals[:,-1], axis=0), np.sqrt(NumEntries))
#ySEs = np.std(BarVals[:,-1], axis=0)
ax2.errorbar(np.array([0]), np.mean(BarVals[:,-1], axis=0), 
             yerr=ySEs, color='gray', ecolor='gray', fmt=ErrorBarFormat, 
             markersize=MarkerSize, label='Shuffled outcomes')

# Calculate Standard Errors of the Means statistics to plot extents of errorbars
#ySEs = np.divide(np.std(ShuffledBarVals[:,1:-1], axis=0), np.sqrt(NumEntries))
ySEs = np.std(ShuffledBarVals[:,1:-1], axis=0)
ax2.errorbar(np.arange(1, NumBarsPerGroup-1), np.mean(ShuffledBarVals[:,1:-1], axis=0), 
             yerr=ySEs, color='lightgray', ecolor='lightgray', fmt=ErrorBarFormat,
             markersize=MarkerSize)

#ax2.errorbar(np.array([0]), ShuffledBarVals[0], 
#             yerr=ShuffledySEs[0], color='gray', ecolor='gray', fmt='.-', 
#             markersize=MarkerSize)

# Set plot axes limits, ticknames, legend location, etc...
ax2.set_xlabel('Units Used for Decoding')
ax2.set_ylabel('Decoding Accuracy (%)')
ax2.set_ylim([0.45, 1.])
ax2.set_yticks(np.arange(.5, .95, .05))
ax2.set_yticklabels(['50','','60','','70','','80','','90'])
ax2.legend(loc='upper right')

# Specify plot labels for x-axis and set accordingly.
xAxisNames = ['All', '50', '40', '30', '20', '10', '1']
ax2.set_xticks(range(0, NumBarsPerGroup, 1))
ax2.set_xticklabels(xAxisNames, ha='center')

# Remove boundary box
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Set title of plot.
ax2.set_title('PLS-DA: 400ms window, slid over 100ms increments')
fig2.suptitle('Trace decoding accuracy averaged over '+ str(NumSessions) +' sessions')
#ax.set_xticklabels(SessionNames, rotation=15, ha="right", va="center")

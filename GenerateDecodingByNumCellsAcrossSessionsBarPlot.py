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

# Integer width to position clusters.
ClusterIncrement = 1.

# The proportion of space to use to separate bar groups.
GroupSepWidthProportion = 0.2

mpl.rc('image', cmap='nipy_spectral')

#PerfMeans = np.empty((NumSamplingDicts,))
#PerfSEs = np.empty((NumSamplingDicts,))
#ShuffledPerfMeans = np.empty((NumSamplingDicts,))
#ShuffledPerfSEs = np.empty((NumSamplingDicts,))

AgrPerfMeans = np.array([])
AgrPerfSEs = np.array([])
AgrShuffledPerfMeans = np.array([])
AgrShuffledPerfSEs= np.array([])

SessionNames = np.array([])

EntriesFilt = np.array([1, 3, 5, 7, 9, -1])
EntriesFilt = EntriesFilt[::-1]

GetAnotherSession = True

while (GetAnotherSession == True):
    
    # Acquire path of workspace to load.
    root = tk.Tk()
    RestoreFilePath = askopenfilename()
    root.withdraw()
    
    # Open workspace.    
    exec(open('./RestoreShelvedWorkspaceScript.py').read())

    # Determine parent directory and filename from complete path.
    drive, path_and_file = os.path.splitdrive(RestoreFilePath)
    path, file = os.path.split(path_and_file)
    
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
        
    GetAnotherSession = tk.messagebox.askyesno(message='Include another session?')
    
    root.withdraw()

BarNames = np.hstack([np.array(np.array(X[EntriesFilt], dtype=int), dtype=str), 'Shuffled'])
BarNames[0] = 'All'
BarVals = np.hstack([AgrPerfMeans, np.array([AgrShuffledPerfMeans[:, 0]]).transpose()])
SEVals = np.hstack([AgrPerfSEs, np.array([AgrShuffledPerfSEs[:, 0]]).transpose()])

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

ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
ax.set_xlabel('Session')
ax.set_ylabel('Decoding Accuracy')
ax.set_ylim([0., 1.])
NiceSessionNames = SessionNames
for i in np.arange(0, SessionNames.shape[0]):
    NiceSessionNames[i] = SessionNames[i][0:10]
    
ax.set_xticks(range(1, GroupPositions.shape[0]+1, 1))
ax.set_xticklabels(NiceSessionNames, rotation=45, ha='center')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_title('Average trace decoding accuracy by no. of cells included across sessions')


fig2, ax2 = plt.subplots(nrows=1, ncols=1)

ShuffledVal = np.mean(BarVals[:,-1])
ax2.errorbar(np.arange(1, NumBarsPerGroup-1), np.mean(BarVals[:,1:-1], axis=0), yerr=np.std(BarVals[:,1:-1], axis=0), color='blue', ecolor='blue', fmt='.-', label='Decoding with subsampled units')
ax2.errorbar(np.array([0]), np.mean(BarVals[:,0], axis=0), yerr=np.std(BarVals[:,0], axis=0), color='black', ecolor='black', fmt='.-', label='Decoding with all units')
ax2.plot(np.array([0, NumBarsPerGroup-2]), np.array([ShuffledVal, ShuffledVal]), 'c--', label='Shuffled performance (all units)')
ax2.set_xlabel('Units Used for Decoding')
ax2.set_ylabel('Decoding Accuracy')
#ax2.set_ylim([0., 1.])
ax2.legend(loc='upper right')
xAxisNames = ['All units', '50 units', '40 units', '30 units', '20 units', '10 units']
ax2.set_xticks(range(0, NumBarsPerGroup, 1))
ax2.set_xticklabels(xAxisNames, ha='center')
ax2.set_title('Trace decoding accuracy by no. of included units averaged across sessions')

#ax.set_xticklabels(SessionNames, rotation=15, ha="right", va="center")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 01:00:07 2019

@author: thugwithyoyo
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

import tkinter as tk
from tkinter.filedialog import askopenfilename

def GetSessionDecodingMaximum():
    
    # Acquire path of workspace to load.
    root = tk.Tk()
    RestoreFilePath = askopenfilename()
    root.withdraw()
    
    # Open workspace.
    #exec(open('./RestoreShelvedWorkspaceScript.py').read())
    try:
        
        exec(open('./RestoreShelvedWorkspaceScript.py').read())
        
    except:
        
        print('Unshelving error.  Will attempt to continue...')

    # Determine parent directory and filename from complete path.
    drive, path_and_file = os.path.splitdrive(RestoreFilePath)
    path, file = os.path.split(path_and_file)
    
    PerfIndices = np.arange(1, Performance.shape[0])
    MaxPerfVal = Performance[0]['performance']
    MaxPerfIndex = 0
    
    SessionName = file[0:19]
            
    for i in PerfIndices:
        
        if Performance[i]['performance'] > MaxPerfVal:
            
            MaxPerfIndex = i
            MaxPerfVal = Performance[i]['performance']
            
    return {
            'MaxPerfIndex': np.array([MaxPerfIndex]),
            'MaxPerfVal': np.array([MaxPerfVal]),
            'MaxPerfErrorBar': np.array([ConfInts[MaxPerfIndex]['performance_SE']]),
            'ShuffledPerfMedian': np.array([EventsShuffled[MaxPerfIndex]['performance_median']]),
            'ShuffledPerfErrorBar': np.array([EventsShuffled[MaxPerfIndex]['performance_SE']]),
            'SessionName': np.array([SessionName])
            }

# Integer width to position clusters.
ClusterIncrement = 1.

# The proportion of space to use to separate bar groups.
GroupSepWidthProportion = 0.25

# Calculate barwidth
BarWidth = ClusterIncrement*(1. - GroupSepWidthProportion)/2.

# Set colorscheme.  Is this necessary?
mpl.rc('image', cmap='nipy_spectral')

SessionPerfMeans = np.empty((0,))
SessionPerfSEs = np.empty((0,))
SessionShuffledPerfMeans = np.empty((0,))
SessionShuffledPerfSEs = np.empty((0,))

SessionNames = np.array([])

#EntriesFilt = np.array([1, 3, 5, 7, 9, -1])
#EntriesFilt = EntriesFilt[::-1]

GetAnotherSession = True

while (GetAnotherSession == True):
    
    SessionDict = GetSessionDecodingMaximum()
        
#    SessionPerfMeans = np.hstack([SessionPerfMeans, np.array([MaxPerfVal])])
#    SessionPerfSEs = np.hstack([SessionPerfSEs, np.array([ConfInts[i]['performance_SE']])])
#    SessionShuffledPerfMeans = np.hstack([SessionShuffledPerfMeans, np.array([EventsShuffled[i]['performance_median']])])
#    SessionShuffledPerfSEs = np.hstack([SessionShuffledPerfSEs, np.array([EventsShuffled[i]['performance_SE']])])
#    
#    SessionNames = np.hstack([SessionNames, SessionName])
    
    SessionPerfMeans = np.hstack([SessionPerfMeans, SessionDict['MaxPerfVal']])
    SessionPerfSEs = np.hstack([SessionPerfSEs, SessionDict['MaxPerfErrorBar']])
    SessionShuffledPerfMeans = np.hstack([SessionShuffledPerfMeans, SessionDict['ShuffledPerfMedian']])
    SessionShuffledPerfSEs = np.hstack([SessionShuffledPerfSEs, SessionDict['ShuffledPerfErrorBar']])
    
    SessionNames = np.hstack([SessionNames, SessionDict['SessionName']])    
    
    root = tk.Tk()
    
    GetAnotherSession = tk.messagebox.askyesno(message='Include another session?')
    
    root.withdraw()

#     
# Extract recording IDs from filenames
NiceSessionNames = SessionNames

for i in np.arange(0, SessionNames.shape[0]):
    
    NiceSessionNames[i] = SessionNames[i][0:10]

SortMap = np.argsort(NiceSessionNames)

# Begin plot generation
(NumSessions, ) = SessionPerfMeans.shape
fig1, ax1 = plt.subplots(nrows=1, ncols=1)

fig1.suptitle('Peak trace decoding accuracy across sessions')

xLocs = np.arange(ClusterIncrement, (NumSessions + 1)*ClusterIncrement, 
                  ClusterIncrement)

# Plot bars and errorbars of OBSERVED peak decoder performance across sessions.
ax1.bar(xLocs - BarWidth/2., SessionPerfMeans[SortMap], 
        color='orange', ecolor='orange', label='Observed outcomes', 
        width=BarWidth)

ax1.errorbar(xLocs - BarWidth/2., SessionPerfMeans[SortMap], 
             yerr=SessionPerfSEs[SortMap], 
             color='black', ecolor='black', fmt=',')

# Plot bars and errorbars of SHUFFLED peak decoder performance across sessions.
ax1.bar(xLocs + BarWidth/2., SessionShuffledPerfMeans[SortMap], 
        color='gray', ecolor='gray', label='Shuffled outcomes', 
        width=BarWidth)

ax1.errorbar(xLocs + BarWidth/2., SessionShuffledPerfMeans[SortMap], 
             yerr=SessionShuffledPerfSEs[SortMap], 
             color='black', ecolor='black', fmt=',')

ax1.set_xlabel('Session')
ax1.set_ylabel('Peak Decoding Accuracy (%)')
ax1.set_ylim([0., 1.])
ax1.set_yticks([0., 0.5, 1.0])
ax1.set_yticklabels(['0', '50', '100'])

    
#ax1.set_xticks(range(1, NumSessions+1, 1))
ax1.set_xticks(xLocs)
#ax1.set_xticklabels(NiceSessionNames[SortMap], rotation=45, ha='center')
ax1.set_xticklabels(np.arange(1,12,1), rotation=0)
ax1.legend(loc='lower right')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.set_title('PLS-DA: 400ms window, slid over 100ms increments')

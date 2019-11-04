#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 01:00:07 2019

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

SessionPerfMeans = np.empty((0,))
SessionPerfSEs = np.empty((0,))
SessionShuffledPerfMeans = np.empty((0,))
SessionShuffledPerfSEs = np.empty((0,))

SessionNames = np.array([])

#EntriesFilt = np.array([1, 3, 5, 7, 9, -1])
#EntriesFilt = EntriesFilt[::-1]

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
    
    PerfIndices = np.arange(1, Performance.shape[0])
    MaxPerfVal = Performance[0]['performance']
    MaxPerfIndex = 0
    
            
    for i in PerfIndices:
        
        if Performance[i]['performance'] > MaxPerfVal:
            
            MaxPerfIndex = i
            MaxPerfVal = Performance[i]['performance']
        
    SessionPerfMeans = np.hstack([SessionPerfMeans, np.array([MaxPerfVal])])
    SessionPerfSEs = np.hstack([SessionPerfSEs, np.array([ConfInts[i]['performance_SE']])])
    SessionShuffledPerfMeans = np.hstack([SessionShuffledPerfMeans, np.array([EventsShuffled[i]['performance_median']])])
    SessionShuffledPerfSEs = np.hstack([SessionShuffledPerfSEs, np.array([EventsShuffled[i]['performance_SE']])])
    
    SessionNames = np.hstack([SessionNames, file[0:19]])
        
    GetAnotherSession = tk.messagebox.askyesno(message='Include another session?')
    
    root.withdraw()


(NumSessions, ) = SessionPerfMeans.shape
fig1, ax1 = plt.subplots(nrows=1, ncols=1)

ax1.errorbar(np.arange(1, NumSessions+1), SessionPerfMeans, yerr=SessionPerfSEs, color='blue', ecolor='blue', fmt='.-', label='Observed outcomes')
ax1.errorbar(np.arange(1, NumSessions+1), SessionShuffledPerfMeans, yerr=SessionShuffledPerfSEs, color='cyan', ecolor='cyan', fmt='.-', label='Shuffled outcomes')

ax1.set_xlabel('Session')
ax1.set_ylabel('Decoding Accuracy')
ax1.set_ylim([0.4, 1.])
NiceSessionNames = SessionNames
for i in np.arange(0, SessionNames.shape[0]):
    NiceSessionNames[i] = SessionNames[i][0:10]
    
ax1.set_xticks(range(1, NumSessions+1, 1))
ax1.set_xticklabels(NiceSessionNames, rotation=45, ha='center')
ax1.legend()
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax1.set_title('Peak trace decoding accuracy across sessions')
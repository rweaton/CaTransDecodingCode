#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 23:40:34 2019

@author: thugwithyoyo
"""

import matplotlib.pyplot as plt
from PeriEventTraceFuncLib import *
import numpy as np

#RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-10-11-37-56/SW_400ms_100ms/2018-12-10-11-37-56_unique_SlidingWindow.dat'
RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-14-11-01-41/JointHemisphereDecoding/SingleArmDirectionTask/2018-12-14-11-01-41_new_unique_SW_LHem.dat'
exec(open('./RestoreShelvedWorkspaceScript.py').read())

RefEventsList = RefEventsDict['RefEventsList']
nEvents = len(RefEventsList)

PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                                  CellFluorTraces_Frame, 
                                                  RefEventsDict, 
                                                  ParamsDict['BoundaryWindow'])

NumTrials, NumCellsByNumSamples = PeriEventExtractorDict['PEA_Array'].shape

(_, NumColumns) = CellFluorTraces_Frame.shape

NumCells = NumColumns - 1

NumSamples = int(NumCellsByNumSamples/NumCells)

RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], ParamsDict['BoundaryWindow'][1], num=NumSamples, endpoint=False)

NumEvents = len(RefEventsList)

AveragedTracesMatrices = np.empty((NumCells, NumSamples, NumEvents))

fig, axes  = plt.subplots(1, NumEvents)

# generate 2 2d grids for the x & y bounds
xv, yv = np.meshgrid(RelTimeVec, np.arange(1, NumCells + 1, 1))
i = 0

for RefEvent in RefEventsList:
    
    TracesFilt = PeriEventExtractorDict['TrialIndicesByEventDict'][RefEvent]
    
    AveragedTraces = np.mean(PeriEventExtractorDict['PEA_Array'][TracesFilt,:],
                                      axis=0)
    
    AveragedTracesMatrices[:,:,i] = np.reshape(AveragedTraces, (NumCells, NumSamples))
    

    #im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
    
    axes[i].pcolormesh(xv, yv, AveragedTracesMatrices[:,:,i])
    
    axes[i].set_xlabel('time relative to reach (sec.)')
    
    if i == 0:
        axes[0].set_ylabel('Cell identity')
    
    axes[i].title.set_text(RefEvent)
    
    i += 1
    
#  Generate difference plot. 
fig2, axs2 = plt.subplots(1,1)

DiffMat = np.diff(AveragedTracesMatrices, axis=-1)

im2 = axs2.pcolormesh(xv, yv, DiffMat[:,:,0])
fig2.colorbar(im2, ax=axs2)
axs2.set_ylabel('Cell identity')
axs2.set_xlabel('time relative to reach (sec.)')

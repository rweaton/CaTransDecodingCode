#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 22:33:00 2019

@author: thugwithyoyo
"""

# Initialize figure
FigureTitle = 'Dual hand task decoder performance by hemisphere'
fig1, axs1 = plt.subplots()
fig1.suptitle(FigureTitle)

LHemPath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-27/JointHemisphereDecoding/2018-12-27-11-35-13_new_unique_SW_LHem.dat'
RestoreFilePath = LHemPath
exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Plot performance and performance control plots
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'red'}

GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
                     axs1, 'fw_sliding')

RHemPath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-27/JointHemisphereDecoding/2018-12-27-11-35-22_new_unique_SW_RHem.dat'
RestoreFilePath = RHemPath
exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Plot performance and performance control plots
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'blue'}

GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
                     axs1, 'fw_sliding')

LandRHemPath = '/home/thugwithyoyo/CaTransDecoding/Output/2018-12-27/JointHemisphereDecoding/2018-12-27_new_unique_SW_LandRHem.dat'
RestoreFilePath = LandRHemPath
exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Plot performance and performance control plots
PlotSpecDict = {'measure': 'performance',
                'measure_median': 'performance_median',
                'measure_CLs': 'performance_CLs',
                'color':'purple'}

GenerateConfIntsPlot(ConfInts, Performance, PlotSpecDict, 
                     axs1, 'fw_sliding')
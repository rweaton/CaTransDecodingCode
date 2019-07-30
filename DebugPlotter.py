#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 16:22:55 2019

@author: thugwithyoyo
"""
fig1, axs1 = plt.subplots()

GenerateConfIntsPlot(ConfInts, Performance, PerformancePlotSpecDict, 
                     axs1, 'fw_sliding')

axs1.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs1.set_ybound(lower=0.4, upper=1.)

# Plot mutual information dependence on increasing peri-event window span
fig2, axs2 = plt.subplots()
#fig2.suptitle(MutInfoPlotSpecDict['measure'])
GenerateConfIntsPlot(ConfInts, Performance, MutInfoPlotSpecDict, 
                     axs2, 'fw_sliding')
axs2.set_xbound(lower=BoundaryWindow[0], upper=BoundaryWindow[1])
axs2.set_ybound(lower=0., upper=1.)

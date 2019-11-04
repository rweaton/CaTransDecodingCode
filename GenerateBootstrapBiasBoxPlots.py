#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:49:20 2019

@author: thugwithyoyo
"""
import matplotlib.pyplot as plt
import numpy as np

TicksToSkip = 4

fig, axs = plt.subplots(ncols=1,nrows=1)
(NumDomains,) = ConfInts.shape
(NumSamples, ) = ConfInts[0]['performance_SortedSimDist'].shape


BiasDistros = np.empty((NumDomains, NumSamples))

positions = np.arange(1, NumDomains+1)

for i in np.arange(0, NumDomains):
    
    BiasDistros[i, :] = (ConfInts[i]['performance_SortedSimDist'] - 
                         Performance[i]['performance'])    

#plt.boxplot(BiasDistros.transpose(), positions=positions, notch=True)

plt.boxplot(BiasDistros.transpose(), notch=True)

xTicks = np.arange(1, NumDomains+1, TicksToSkip)

axs.set_xticks(xTicks)
axs.set_xticklabels(np.array(xTicks, dtype=str))

axs.set_ylabel(r'Bootstrapped perf. - observed perf.: t($\hat F$*) - t{$\hat F$}')
axs.set_xlabel('Sliding window no.')
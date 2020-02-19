#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:16:54 2019

@author: thugwithyoyo
"""

import numpy as np
import PartialLeastSquares as pls

### Set constants ###
epsilon = 3.
nChoices = 150
t_window = [-1., 1.]
nSamples = 100

AmpVec = np.array([1., 1., 1.])
#AmpVec = np.array([1., 2., 3.])

#FreqVec = np.array([1., 1.5, 2.])
FreqVec = np.array([1., 1., 1.])

PhaseVec = np.array([0., np.pi/4., np.pi/2.])
#PhaseVec = np.array([0., 0., 0.])

# Set parameters for PLS
nLatents = 2

ChoiceVec = np.round(2*np.random.rand(nChoices))
t_vec = np.linspace(t_window[0], t_window[1], num=nSamples, endpoint=True)

# Initialize array to contain data
X = np.empty((nChoices, t_vec.size))

# Generate frequency-varied seed arrays
ZeroSeed = np.sin(2*np.pi**t_vec)
OneSeed = np.sin(2*np.pi*FreqVec[1]*t_vec)
TwoSeed = np.sin(2*np.pi*FreqVec[2]*t_vec)

# Generate frequency-, phase- and amplitude-varied seed arrays
ZeroSeed = AmpVec[0]*np.sin(2.*np.pi*FreqVec[0]*t_vec + PhaseVec[0])
OneSeed = AmpVec[1]*np.sin(2.*np.pi*FreqVec[1]*t_vec + PhaseVec[1])
TwoSeed = AmpVec[2]*np.sin(2.*np.pi*FreqVec[2]*t_vec + PhaseVec[2])
       
for i in np.arange(0, ChoiceVec.size):
    
   if (ChoiceVec[i] == 0):
       
       X[i, :] = ZeroSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)
       
   if (ChoiceVec[i] == 1):
       
       X[i, :] = OneSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)
       
   if (ChoiceVec[i] == 2):
       
       X[i, :] = TwoSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)       
       
B, B_0, _ = pls.PLS1(X, ChoiceVec, nLatents)

Predicted = X @ B + B_0

#nIncorrects = np.sum(np.abs(np.squeeze(np.round(Predicted)) - ChoiceVec))
nIncorrects = np.sum(np.squeeze(np.round(Predicted)) != ChoiceVec)
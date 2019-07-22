#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 13:16:54 2019

@author: thugwithyoyo
"""

import numpy as np
import PartialLeastSquares as pls

### Set constants ###
epsilon = 5.
nChoices = 150
t_window = [-1., 1.]
nSamples = 100

# Set parameters for PLS
nLatents = 5

ChoiceVec = np.round(2*np.random.rand(nChoices))
t_vec = np.linspace(t_window[0], t_window[1], num=nSamples, endpoint=True)

# Initialize array to contain data
X = np.empty((nChoices, t_vec.size))

for i in np.arange(0, ChoiceVec.size):
    
   ZeroSeed = np.sin(2*np.pi*t_vec)
   OneSeed = np.sin(2*1.5*np.pi*t_vec)
   TwoSeed = np.sin(2*2.0*np.pi*t_vec)
   
   if (ChoiceVec[i] == 0):
       
       X[i, :] = ZeroSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)
       
   if (ChoiceVec[i] == 1):
       
       X[i, :] = OneSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)
       
   if (ChoiceVec[i] == 2):
       
       X[i, :] = TwoSeed + epsilon*(np.random.rand(t_vec.size) - 0.5)       
       
B, B_0 = pls.PLS1(X, ChoiceVec, nLatents)

Predicted = X @ B + B_0

nIncorrects = np.sum(np.abs(np.squeeze(np.round(Predicted)) - ChoiceVec))
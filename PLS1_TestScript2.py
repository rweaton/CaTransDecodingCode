#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 18:39:07 2020

@author: thugwithyoyo
"""

import numpy as np
import PartialLeastSquares as pls
import hwfunclib as HFL
import matplotlib.pyplot as plt

# Set parameters for PLS
nLatents = 3

############################## Set parameters #################################
# Set amplitude of noise to be added to seed predictors
epsilon = 0.5

# Set the total number of predictor-outcome (choice) pairs to be generated
nChoices = 300

# Specify the number of unique predictor types on which to test the classifier
nSeeds = 5

# Specify the type of underlying signal of seed vectors to which noise will be
# added to construct predictors
#SignalType = 'sine'
#SignalType = 'gaussian'
SignalType = 'multishape'

# Set the domain boundaries and sampling rate of the time window corresponding
# to predictor vectors. 
t_window = [-1., 1.]
nTimeSamples = 100

# Specify the fraction of total predictors to be used to train the model.  The
# remaining predictors will be used to test performance of the model.
TrainingFrac = 0.8

############# end set parameters ###########

# Generate the vector of time values to correpond with predictor vectors
t_vec = np.linspace(t_window[0], t_window[1], num=nTimeSamples, endpoint=False)

############ Parameters for sine wave predictor seed functions ################    
# Build array of amplitudes
AmpVec = np.ones((nSeeds,), dtype=float)
#AmpVec = np.arange(1, nSeeds + 1, 1, dtype=float)

# Build array of frequencies
FreqVec = np.ones((nSeeds,), dtype=float)
#FreqVec = np.arange(1, nSeeds + 1, 1)

# Build array of phases
#PhaseVec = np.zeros((nSeeds,), dtype=float)
PhaseVec = (np.pi/4.)*np.arange(0, nSeeds, 1)


############## Parameters for Gaussian predictor seed vectors. ################
# Build a shifting mean vector
CombSpacing = np.floor(nTimeSamples/nSeeds)
Comb = np.arange(np.round(CombSpacing/2.), nTimeSamples, CombSpacing, dtype=int)
MeanVec = t_vec[Comb]

# Generate uniform-valued vector of standard deviations
StdDevVec = (1./3.)*np.ones((nSeeds,), dtype=float)

#################### End parameters for seed functions ########################

# Build array of model target values that are symmetric about zero.
if np.mod(nSeeds, 2.)  == 1.:
    
    TargetValues = np.arange(-np.floor(nSeeds/2.), np.ceil(nSeeds/2.), 1)

if np.mod(nSeeds, 2.) == 0.:
    
    TargetValues = np.arange(-np.floor(nSeeds/2.), np.ceil(nSeeds/2.), 1) + 0.5

# Randomly sample, with replacement, from the list of target types n=nChoices
# times to generate the list of desired outcomes to be used to train the model.
# Each outcome will correspond to an individual seed vector.
ChoiceVec = np.random.choice(TargetValues, size=(nChoices,))

# Initialize array to contain the seed vectors
SeedVecs = np.empty((nSeeds, t_vec.size))

# Assemble array of seed vectors
for i in np.arange(0, nSeeds):
    
    if SignalType == 'sine':
        
        SeedVecs[i,:] = AmpVec[i]*np.sin(2.*np.pi*FreqVec[i]*t_vec + PhaseVec[i])
    
    if SignalType == 'gaussian':
        
        SeedVecs[i, :] = AmpVec[i]*np.exp(-(1./2.)*((t_vec - MeanVec[i])/(StdDevVec[i]))**2)

    if SignalType == 'multishape':
        
        FourierSeries = np.zeros_like(t_vec)
        
        for n in np.arange(1, 2*i + 2,  2):
            
            FourierSeries = FourierSeries + (4./np.pi)*(1./n)*np.sin(2*n*np.pi*FreqVec[i]*t_vec)
            
        SeedVecs[i, :] = AmpVec[i]*FourierSeries

# Initialize array to contain predictor vectors
X = np.empty((nChoices, t_vec.size))
    
# Generate list of indexes to locate target values
NumTargetValues = TargetValues.size
TargValIndices = np.arange(0, NumTargetValues)
    
# Populate the array of predictor vectors; each is the seed vector, with added 
# noise, that corresponds to the particular choice for each iteration.
for i in np.arange(0, ChoiceVec.size):
    
   j = TargValIndices[TargetValues == ChoiceVec[i]]
       
   X[i, :] = SeedVecs[j, :] + epsilon*(np.random.rand(t_vec.size) - 0.5)
   
# Construct slice definitions to extract training and testing sets of predictors
TrainingLimitIndex = int(np.round(TrainingFrac*nChoices))
TrainingSubset = slice(0, TrainingLimitIndex, 1)

# Run partial least squares classification using the training subset
B, B_0, _ = pls.PLS1(X[TrainingSubset, :], ChoiceVec[TrainingSubset], nLatents)

# Extract the testing subset of predictors. Run testing subset in trained model.
TestingSubset = slice(TrainingLimitIndex + 1, -1)
Predicted = X[TestingSubset, :] @ B + B_0

# Calculate Euclidean distance
EuclideanDiffVals = np.empty_like(TargetValues)

# Initialize the confusion matrix
ConfusionMatrix = np.zeros((NumTargetValues, NumTargetValues))

# Initialize the array to contain model predictions
PredictedTargets = np.empty_like(ChoiceVec[TestingSubset])

# Initialize the array to contain actual corresponding targets
TrueTargets = np.empty_like(ChoiceVec[TestingSubset])

# For each trail of the test set predict target type by selecting the closest
# value to model output.  Compare with correct target list.  Construct
# confusion matrix.
for i in np.arange(0, ChoiceVec[TestingSubset].size):
    
  for k in TargValIndices:
        
    EuclideanDiffVals[k] = HFL.EuclideanDistanceCalculator(Predicted[i], 
                     TargetValues[k])

  # Acquire index of closest event.  This is the model's prediction.
  PredictedEventIndex = np.argmin(EuclideanDiffVals)
  
  # Obtain the value of the actual target.
  TrueEventIndex = TargValIndices[TargetValues == ChoiceVec[TestingSubset][i]][0]
        
  # Increment the cooresponding element of the confusion matrix (e.g. correct
  # prediction (on diagonal), or incorrect prediction (off diagonal))
  ConfusionMatrix[TrueEventIndex, PredictedEventIndex] += 1        
        
  # Record the prediction.
  PredictedTargets[i] = TargetValues[PredictedEventIndex]
  
  # Report the actual corresponding target for the current predictor
  TrueTargets[i] = TargetValues[TrueEventIndex]
  
# Plot output predictor vectors in subplots of each target type.  Grey traces
# show predictors that were predicted correctly.  Red show predictors that
# were supposed to be classified in that set, but were misclassified by 
# the model.
fig, axs = plt.subplots(nrows=nSeeds, ncols=2)

# Iterate through list of target types.  For each type, plot correctly-classified
# traces in gray.  Plot incorrectly classified traces in red.
for k in TargValIndices:
  
    PredictedTargValFilt = (PredictedTargets == TargetValues[k])
    TrueTargetsFilt = (TrueTargets == TargetValues[k])
    
    CorrectFilt = (TrueTargetsFilt & PredictedTargValFilt)
  
    IncorrectFilt = (TrueTargetsFilt & ~PredictedTargValFilt)
    
    if np.sum(CorrectFilt) > 0:
        
        axs[k,0].plot(t_vec, X[TestingSubset, :][CorrectFilt, :].transpose(), 
           color='gray')
  
    if np.sum(IncorrectFilt) > 0:
      
        axs[k,1].plot(t_vec, X[TestingSubset, :][IncorrectFilt, :].transpose(), 
           color='red')
  
    if k != TargValIndices[-1]:

        axs[k,0].set_xticklabels([])
        axs[k,1].set_xticklabels([])

#nIncorrects = np.sum(np.abs(np.squeeze(np.round(Predicted)) - ChoiceVec))
#nIncorrects = np.sum(np.squeeze(np.round(Predicted)) != ChoiceVec[TestingSubset])
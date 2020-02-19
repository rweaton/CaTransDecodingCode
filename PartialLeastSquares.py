#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 19 08:48:07 2019

@author: thugwithyoyo
"""
import numpy as np

def PLS1(X, y, l):
    
    # Partial Least Squares Regression
    # PLS1 pseudocode from Wikipedia
    # 
    #### Model ####
    #
    # The general underlying model of multivariate PLS is
    #
    #   X = T P.transpose + E
    #   Y = U Q.transpose + F
    #
    # where X is an n x m matrix of predictors, y is an n x p matrix of 
    # responses; T and U are n x l matrices that are, respectively, projections
    # of X (the X score, component or factor matrix) and projections of Y (the
    # Y scores); P and Q are, respectively, m x l and p x l orthogonal loading
    # matrices; and matrices E and F are the error terms, assumed to be
    # independent and indentically distributed random normal variables.  The
    # decompositions of X and Y are made so as to maximize the covariance
    # between T and U.
    #
    #### Algorithm ####
    #
    # This routine constructs a linear regression estimate between X and Y as
    #
    #   Y = X B + B_0.
    #   
    #
    #### Argument definitions ####
    #
    # X: (n x m) data matrix-- or more accurately called the matrix of 
    # predictors.  Value n is the total number of trials, m is the
    # product of no. sources (e.g. cells) with the no. of peri-event samples (e.g.
    # snippets of calcium fluorescence level over time relative to a behavioral
    # event of some kind).
    #
    # y: (n x p) response vector (in this case p = 1).  Value n is the total 
    # number of trials.  A response to each trial is either scaler or binary 
    # valued.
    # 
    # l: the number of latent sources to estimate.
    #
    #### Begin routine ####
    #
    # Convert all values to numpy arrays
    X = np.array(X)
    
    # If y is a  1-D array (likely) convert to a column vector array.
    if (len(y.shape) <  2):
        
        y = np.array([y]).transpose()  # make y a column vector
    
    # Determine size of data matrix X to initialize weights, score and loadings
    # matrices
    (nTrials, nSourcesBySamples) = X.shape
    
    # Initialize weight and loading matrices as well as score and "mixing" arrays
    W = np.empty((nSourcesBySamples, l))

    T_Scores = np.empty((nTrials, l))

    q_Loadings = np.empty((1, l))

    P_Loadings = np.empty((nSourcesBySamples, l))
#    M = np.empty((nTrials, l))
    
    # Initialize an array of indices for iterating across latent sources
    l_indices = np.arange(0, l)
       
    # Set pre-loop variables for iteration number 0
    X_Current = X
    
    # Compute first weight vector: the normalized inner product of the first 
    # vector of data matrix with the output vector.  Recall that "@" is numpy's
    # matrix multiplication operator.
    W[:, 0] = np.squeeze(X.transpose() @ y / np.linalg.norm(X.transpose() @ y))
    
    # Iterate to calculate latent sources
    for k in l_indices:
        
        # Take matrix product of current data vector and weight vector to
        # compute current score vector.
        T_Scores[:, k] = X_Current @ W[:, k]
        
        # Compute squared magnitude of T_scores vector for subsequent scaling
        T_Scores_MagSquared = T_Scores[:, k].transpose() @ T_Scores[:, k]
        
        # Scale the current score vector
        T_Scores[:, k] = T_Scores[:, k] / T_Scores_MagSquared
        
        # Take matrix product of current data vector with the current scores
        # vector to calculate current loadings vector of the input model 
        # equation.
        P_Loadings[:, k] = X_Current.transpose() @ T_Scores[:, k]
        
        # Take matrix product of output vector with the scores vector to 
        # compute the current loadings vector on the output model equation.
        # Note dimensionality of the resulting array must reduced to write
        # the result into the matrix of  Q of loading vectors.
        q_Loadings[0, k] = np.squeeze(y.transpose() @ T_Scores[:, k])
        
        # Check the for-loop termination condition.  If met, set the number
        # of latents to include, l, to the value of the current iteration
        # index, k, so as to skip further optimization (i.e. learning updates in
        # the next step).
        if (q_Loadings[0, k] == 0.):
            
            l = k  # exit for loop
            
        # If the current q_Loadings are sub-optimal (i.e. q_Loadings[0, k] != 0.)
        # then update the current data vector (via NIPALS iterative learning 
        # rule) and compute the weight ("mixing") vector to be used on the 
        # next iteration of the for-loop.
        if (k < (l - 1)):
            
            X_Current = X_Current - T_Scores_MagSquared * (
                    np.array([T_Scores[:, k]]).transpose() @ np.array([P_Loadings[:, k]]))
            
            W[:, k + 1] = np.squeeze(X_Current.transpose() @ y)
    
    # Compute partial least squares regression estimates B and B_0 for Y = X @ B + B_0.
    # Note that this calculation involves a matrix inversion.         
    B = W @ np.linalg.inv(P_Loadings.transpose() @ W) @ q_Loadings.transpose()
    
    B_0 = q_Loadings[0, 0] - P_Loadings[:, 0].transpose() @ B
    
    ModelDict = {
                  'W': W,
                  'T_Scores': T_Scores,
                  'P_Loadings': P_Loadings,
                  'q_Loadings': q_Loadings
                 }
    
    return B, B_0, ModelDict
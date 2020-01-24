#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 10:32:05 2019

@author: thugwithyoyo
"""
import CalciumImagingFluorProcessing as CIFP

PathToFluorFile = '/home/thugwithyoyo/CaTransDecoding/CalciumImagingData/2018-11-26/2018-11-26-11-45-46_new_unique_C.json'

CentroidFrame = CIFP.CellCentroidsFromJSON(PathToFluorFile)
CellFluorTraces_Frame = CIFP.FrameFromJSON(PathToFluorFile)
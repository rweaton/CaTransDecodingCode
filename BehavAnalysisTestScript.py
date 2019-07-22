#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 22:07:25 2019

@author: thugwithyoyo
"""
import CalciumImagingBehaviorAnalysis as CIBA

PathToFile = '/home/thugwithyoyo/Desktop/MAE298_Project/CalciumImagingData/20181210_both_001_1-01Dat.csv'

BehavDict = CIBA.CinePlexCSV_parser(PathToFile)
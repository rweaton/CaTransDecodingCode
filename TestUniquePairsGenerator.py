#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 13:52:26 2019

@author: thugwithyoyo
"""
import PeriEventTraceFuncLib

IntList = np.arange(0,10)

(NumEntries,) = IntList.shape
 
CellIDs = np.empty((NumEntries), dtype=object)

for j in IntList:
    
    CellIDs[j] = 'C' + str(j).zfill(3)
    
CellIDs2 = np.hstack([CellIDs,'C006','C001'])

UniqueIndexPairsList, UniquePairsList = UniquePairsGenerator(CellIDs2)
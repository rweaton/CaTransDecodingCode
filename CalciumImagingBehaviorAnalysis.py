# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:33:38 2019

@author: Ryan Eaton
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Aquire path of the directory that contains this script.  This is so 
# the routine can navigate back after changing directories if need be.  
ScriptDir = os.getcwd()
PathToFile = ''
CellFluorTraces_Frame = pd.read_csv(PathToFile, header=0)

# Generate an index list for later referencing
Indices = np.arange(0, CellFluorTraces_Frame['#'].size)

# Define a dictionary of key name pairs for referencing columns in the dataframe.
# 
# ### Description of nomenclature.  ###
# Here "M6" and "M7" stand for Marker 6 (right hand) and Marker 7 (left hand)
# respectively.  "T0", "T1" and "CT" stand for Target 0 (animal's right), 
# Target 1 (animal's left) and Center Target respectively.  Tags "in" and "out"
# refer to when the relevant marker is present or absent respectively within
# the relevant target. Finally, "CumLen" stands for cumulative length in the
# present target.
ColumnNames = {'IndexNumber':'#',	
               'Timestamp':'Timestamp',
               'FrameNumber':'Frame_Number', 
               'M6_xCoord':'X1.1_pix',
               'M6_yCoord':'Y1.1_pix',
               'M7_xCoord':'X2.1_pix',	
               'M7_yCoord':'Y2.1_pix',
               'M6T0_in':'EV1.1',
               'M6T0_in_CumLen':'EV1.1_Track_Length_pix',
               'M6T0_out':'EV1.2',
               'M6T0_out_CumLen':'EV1.2_Track_Length_pix',
               'M6T1_in':'EV1.3',
               'M6T1_in_CumLen':'EV1.3_Track_Length_pix',
               'M6T1_out':'EV1.4',
               'M6T1_out_CumLen':'EV1.4_Track_Length_pix',
               'M6CT_in':'EV1.5',
               'M6CT_in_CumLen':'EV1.5_Track_Length_pix',
               'M6CT_out':'EV1.6',
               'M6CT_out_CumLen':'EV1.6_Track_Length_pix',
               'M7T0_in':'EV1.7',
               'M7T0_in_CumLen':'EV1.7_Track_Length_pix',
               'M7T0_out':'EV1.8',
               'M7T0_out_CumLen':'EV1.8_Track_Length_pix',
               'M7T1_in':'EV1.9',
               'M7T1_in_CumLen':'EV1.9_Track_Length_pix',
               'M7T1_out':'EV1.10',
               'M7T1_out_CumLen':'EV1.10_Track_Length_pix',
               'M7CT_in':'EV1.11',
               'M7CT_in_CumLen':'EV1.11_Track_Length_pix',
               'M7CT_out':'EV1.12',
               'M7CT_out_CumLen':'EV1.12_Track_Length_pix'}

# Extract Marker 6 T0 entry events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M6T0_in']])== 1))
M6T0_Entry_ind = Indices[Filt]
M6T0_Entry_ts = CellFluorTraces_Frame[ColumnNames['M6T0_in']][M6T0_Entry_ind]

# Extract Marker 6 T1 entry events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M6T1_in']])== 1))
M6T1_Entry_ind = Indices[Filt]
M6T1_Entry_ts = CellFluorTraces_Frame[ColumnNames['M6T1_in']][M6T1_Entry_ind]

# Extract Marker 6 CT exit events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M6CT_out']])== 1))
M6CT_Exit_ind = Indices[Filt]
M6CT_Exit_ts = CellFluorTraces_Frame[ColumnNames['M6CT_out']][M6CT_Exit_ind]

# Extract Marker 7 T0 entry events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M7T0_in']])== 1))
M7T0_Entry_ind = Indices[Filt]
M7T0_Entry_ts = CellFluorTraces_Frame[ColumnNames['M7T0_in']][M7T0_Entry_ind]

# Extract Marker 7 T1 entry events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M7T1_in']])== 1))
M7T1_Entry_ind = Indices[Filt]
M7T1_Entry_ts = CellFluorTraces_Frame[ColumnNames['M7T1_in']][M6T1_Entry_ind]

# Extract Marker 7 CT exit events
Filt = np.hstack((False, np.diff(CellFluorTraces_Frame[ColumnNames['M7CT_out']])== 1))
M7CT_Exit_ind = Indices[Filt]
M7CT_Exit_ts = CellFluorTraces_Frame[ColumnNames['M7CT_out']][M7CT_Exit_ind]

# Write a function that tests if each entry t_2 in tlist2 is within the range 
# [t_1 + tol_lo, t_1 + tol_hi) for each entry in tlist1.  Output is an object array of size 
# tlist1.size having each field contain a list of indices that point to all entries
# in tlist2 that lie within the range. Note, for t + tol_lo to precede t, tol_lo must
# be negative-valued.
def EventComparator(tlist1, tlist2, tol_window):
  
  # Extract low tolerance and high tolerance boundaries from the tol_window tuple.
  (tol_lo, tol_hi) = tol_window
  
  # Generate index lists for subsequent iteration and filtering procedures.
  tlist2_ind = np.arange(0, tlist2.size)
  tlist1_ind = np.arange(0, tlist1.size)
  
  # Pre-allocate array dimensions for output array.
  output = np.nan*np.ones(tlist1.size,1)
  
  for i in tlist1_ind:
    
    filt = (tlist2 >= tlist1[i] + tol_lo) and (tlist2 < tlist1[i] + tol_hi)
    output[i] = np.array([tlist2[filt], tlist2_ind[filt]])
   
  return output
  

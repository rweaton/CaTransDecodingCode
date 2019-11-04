# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 13:33:38 2019

@author: Ryan Eaton
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import collections

def CinePlexFromMatJSON_parser(PathToBehavFile):
    
    ScriptDir = os.getcwd()
   
    with open(PathToBehavFile) as  f:
        BehavDataDict = json.load(f)


    BehaviorTraces_Frame = pd.DataFrame.from_dict(BehavDataDict, 
                                                  orient='columns')
    # Generate an index list for later referencing
    Indices = np.arange(0, BehaviorTraces_Frame.index.size)
    
    # Define a dictionary of key name pairs for referencing columns in the dataframe.
    # 
    # ### Description of nomenclature.  ###
    # Here "M6" and "M7" stand for Marker 6 (right hand) and Marker 7 (left hand)
    # respectively.  "T0", "T1" and "CT" stand for Target 0 (animal's right), 
    # Target 1 (animal's left) and Center Target respectively.  Tags "in" and "out"
    # refer to when the relevant marker is present or absent respectively within
    # the relevant target. Finally, "CumLen" stands for cumulative length in the
    # present target.
#    ColumnNames = {	
#                   'Timestamp':'Timestamp',
#                   'FrameNumber':'Frame_Number', 
#                   'M6_xCoord':'X1_1_pix',
#                   'M6_yCoord':'Y1_1_pix',
#                   'M7_xCoord':'X2_1_pix',	
#                   'M7_yCoord':'Y2_1_pix',
#                   'M6T0_in':'EV1_1',
#                   'M6T0_in_CumLen':'EV1_1_Track_Length_pix',
#                   'M6T0_out':'EV1_2',
#                   'M6T0_out_CumLen':'EV1_2_Track_Length_pix',
#                   'M6T1_in':'EV1_3',
#                   'M6T1_in_CumLen':'EV1_3_Track_Length_pix',
#                   'M6T1_out':'EV1_4',
#                   'M6T1_out_CumLen':'EV1_4_Track_Length_pix',
#                   'M6CT_in':'EV1_5',
#                   'M6CT_in_CumLen':'EV1_5_Track_Length_pix',
#                   'M6CT_out':'EV1_6',
#                   'M6CT_out_CumLen':'EV1_6_Track_Length_pix',
#                   'M7T0_in':'EV1_7',
#                   'M7T0_in_CumLen':'EV1_7_Track_Length_pix',
#                   'M7T0_out':'EV1_8',
#                   'M7T0_out_CumLen':'EV1_8_Track_Length_pix',
#                   'M7T1_in':'EV1_9',
#                   'M7T1_in_CumLen':'EV1_9_Track_Length_pix',
#                   'M7T1_out':'EV1_10',
#                   'M7T1_out_CumLen':'EV1_10_Track_Length_pix',
#                   'M7CT_in':'EV1_11',
#                   'M7CT_in_CumLen':'EV1_11_Track_Length_pix',
#                   'M7CT_out':'EV1_12',
#                   'M7CT_out_CumLen':'EV1_12_Track_Length_pix'
#                   }
    # Samantha's New Column name style
    ColumnNames = {	
                   'Timestamp':'Timestamp',
                   'FrameNumber':'Frame_number', 
                   'M6_xCoord':'RH_x_pos',
                   'M6_yCoord':'RH_y_pos',
                   'M7_xCoord':'LH_x_pos',	
                   'M7_yCoord':'LH_y_pos',
                   'M6T0_in':'RH_zone1',
                   'M6T0_in_CumLen':'EV1_1_Track_Length_pix',
                   'M6T0_out':'EV1_2',
                   'M6T0_out_CumLen':'EV1_2_Track_Length_pix',
                   'M6T1_in':'RH_zone2',
                   'M6T1_in_CumLen':'EV1_3_Track_Length_pix',
                   'M6T1_out':'EV1_4',
                   'M6T1_out_CumLen':'EV1_4_Track_Length_pix',
                   'M6CT_in':'RH_homezone',
                   'M6CT_in_CumLen':'EV1_5_Track_Length_pix',
                   'M6CT_out':'EV1_6',
                   'M6CT_out_CumLen':'EV1_6_Track_Length_pix',
                   'M7T0_in':'LH_zone1',
                   'M7T0_in_CumLen':'EV1_7_Track_Length_pix',
                   'M7T0_out':'EV1_8',
                   'M7T0_out_CumLen':'EV1_8_Track_Length_pix',
                   'M7T1_in':'LH_zone2',
                   'M7T1_in_CumLen':'EV1_9_Track_Length_pix',
                   'M7T1_out':'EV1_10',
                   'M7T1_out_CumLen':'EV1_10_Track_Length_pix',
                   'M7CT_in':'LH_homezone',
                   'M7CT_in_CumLen':'EV1_11_Track_Length_pix',
                   'M7CT_out':'EV1_12',
                   'M7CT_out_CumLen':'EV1_12_Track_Length_pix'
                   }
    
    # Assuming that the imported data has already been corrected for
    # pause delay.
    
    # Extract Marker 6 T0 entry events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T0_in']])== 1))
    M6T0_Entry_ind = Indices[Filt]
    M6T0_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T0_Entry_ind]
        
    # Extract Marker 6 T0 exit events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T0_in']])== -1))
    M6T0_Exit_ind = Indices[Filt]
    M6T0_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T0_Exit_ind]
    
    # Extract Marker 6 T1 entry events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T1_in']])== 1))
    M6T1_Entry_ind = Indices[Filt]
    M6T1_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T1_Entry_ind]
    
    # Extract Marker 6 T1 exit events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T1_in']])== -1))
    M6T1_Exit_ind = Indices[Filt]
    M6T1_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T1_Exit_ind]    
    
    # Extract Marker 6 CT exit events
    #Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6CT_out']])== 1))
    #M6CT_Exit_ind = Indices[Filt]
    #M6CT_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6CT_Exit_ind]

    # Extract Marker 7 T0 entry events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T0_in']])== 1))
    M7T0_Entry_ind = Indices[Filt]
    M7T0_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T0_Entry_ind]
    
    # Extract Marker 7 T0 exit events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T0_in']])== -1))
    M7T0_Exit_ind = Indices[Filt]
    M7T0_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T0_Exit_ind]    
    
    # Extract Marker 7 T1 entry events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T1_in']])== 1))
    M7T1_Entry_ind = Indices[Filt]
    M7T1_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T1_Entry_ind]
    
    # Extract Marker 7 T1 entry events
    Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T1_in']])== -1))
    M7T1_Exit_ind = Indices[Filt]
    M7T1_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T1_Exit_ind]
    
    # Extract Marker 7 CT exit events
    #Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7CT_out']])== 1))
    #M7CT_Exit_ind = Indices[Filt]
    #M7CT_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7CT_Exit_ind]
 
    return {
            'M6T0_Entry_ind' : M6T0_Entry_ind,
            'M6T0_Entry_ts': M6T0_Entry_ts,
            'M6T0_Exit_ind' : M6T0_Exit_ind,
            'M6T0_Exit_ts' : M6T0_Exit_ts,
            'M6T1_Entry_ind' : M6T1_Entry_ind,
            'M6T1_Entry_ts' : M6T1_Entry_ts,
            'M6T1_Exit_ind' : M6T1_Exit_ind,
            'M6T1_Exit_ts'  : M6T1_Exit_ts,
            #'M6CT_Exit_ind' : M6CT_Exit_ind,
            #'M6CT_Exit_ts' : M6CT_Exit_ts,
            'M7T0_Entry_ind' : M7T0_Entry_ind,
            'M7T0_Entry_ts' : M7T0_Entry_ts,
            'M7T0_Exit_ind' : M7T0_Exit_ind,            
            'M7T0_Exit_ts' : M7T0_Exit_ts,           
            'M7T1_Entry_ind' : M7T1_Entry_ind,
            'M7T1_Entry_ts' : M7T1_Entry_ts,
            'M7T1_Exit_ind' : M7T1_Exit_ind,
            'M7T1_Exit_ts'  : M7T1_Exit_ts,            
            #'M7CT_Exit_ind' : M7CT_Exit_ind,
            #'M7CT_Exit_ts' : M7CT_Exit_ts
           }
 
def CinePlexCSV_parser(PathToFile):
  
  # Function takes a CSV file of present/absent events, converts changes
  # of these states into timestamps, and returns each of the relevant
  # lists of entry and exit events in a dictionary for subsequent
  # processing.
  
  # Aquire path of the directory that contains this script.  This is so 
  # the routine can navigate back after changing directories if need be.  
  ScriptDir = os.getcwd()
  #PathToFile = ''
  
  BehaviorTraces_Frame = pd.read_csv(PathToFile, header=0)
  
#  try:
#      
#      BehaviorTraces_Frame = pd.read_csv(PathToFile, header=0)
#      
#  except:
#      
#      print('File: ' + PathToFile + 'generated a pd.read_csv error...\n')
#      
#  return BehaviorTraces_Frame

  # Generate an index list for later referencing
  Indices = np.arange(0, BehaviorTraces_Frame['#'].size)

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
  
  # Identify "pre-pause" delay in timestamp record
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['Timestamp']].values) >
                    1.1*(BehaviorTraces_Frame[ColumnNames['Timestamp']].values[-1] -
                     BehaviorTraces_Frame[ColumnNames['Timestamp']].values[-2])))
  
  if np.sum(Filt) > 0:
      
      StartIndex = Indices[Filt][0]
      
  elif np.sum(Filt) == 0:
      
      StartIndex = 0
      
  # Subtract out pre-pause delay on startup
  BehaviorTraces_Frame[ColumnNames['Timestamp']] = \
      BehaviorTraces_Frame[ColumnNames['Timestamp']].values - \
      BehaviorTraces_Frame[ColumnNames['Timestamp']].values[StartIndex]
  
  # Remove rows from dataframe that correspond to the in-pause time.
  RowsToKeepFilt = (BehaviorTraces_Frame[ColumnNames['Timestamp']].values >= 0)
  BehaviorTraces_Frame = BehaviorTraces_Frame.iloc[RowsToKeepFilt,:]
  
  # Generate a list of indices for the resulting filtered table.
  Indices = np.arange(0, BehaviorTraces_Frame['#'].size)
  
  # Reset indices of dataframe
  BehaviorTraces_Frame.index = Indices
                                  
  # Extract Marker 6 T0 entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T0_in']])== 1))
  M6T0_Entry_ind = Indices[Filt]
  M6T0_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T0_Entry_ind]

  # Extract Marker 6 T0 exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T0_out']])== 1))
  M6T0_Exit_ind = Indices[Filt]
  M6T0_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T0_Exit_ind]
  
  # Extract Marker 6 T1 entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T1_in']])== 1))
  M6T1_Entry_ind = Indices[Filt]
  M6T1_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T1_Entry_ind]
  
  # Extract Marker 6 T1 exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6T1_out']])== 1))
  M6T1_Exit_ind = Indices[Filt]
  M6T1_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6T1_Exit_ind]  

  # Extract Marker 6 CT entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6CT_in']])== 1))
  M6CT_Entry_ind = Indices[Filt]
  M6CT_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6CT_Entry_ind]
  
  # Extract Marker 6 CT exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M6CT_out']])== 1))
  M6CT_Exit_ind = Indices[Filt]
  M6CT_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M6CT_Exit_ind]

  # Extract Marker 7 T0 entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T0_in']])== 1))
  M7T0_Entry_ind = Indices[Filt]
  M7T0_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T0_Entry_ind]

  # Extract Marker 7 T0 exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T0_out']])== 1))
  M7T0_Exit_ind = Indices[Filt]
  M7T0_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T0_Exit_ind]
  
  # Extract Marker 7 T1 entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T1_in']])== 1))
  M7T1_Entry_ind = Indices[Filt]
  M7T1_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T1_Entry_ind]

  # Extract Marker 7 T1 exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7T1_out']])== 1))
  M7T1_Exit_ind = Indices[Filt]
  M7T1_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7T1_Exit_ind]
 
  # Extract Marker 7 CT entry events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7CT_in']])== 1))
  M7CT_Entry_ind = Indices[Filt]
  M7CT_Entry_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7CT_Entry_ind]
  
  # Extract Marker 7 CT exit events
  Filt = np.hstack((False, np.diff(BehaviorTraces_Frame[ColumnNames['M7CT_out']])== 1))
  M7CT_Exit_ind = Indices[Filt]
  M7CT_Exit_ts = BehaviorTraces_Frame[ColumnNames['Timestamp']][M7CT_Exit_ind]
  
  return {# Right hand enters right target
          'M6T0_Entry_ind' : M6T0_Entry_ind,
          'M6T0_Entry_ts': M6T0_Entry_ts,
          # Right hand exits right target
          'M6T0_Exit_ind' : M6T0_Exit_ind,
          'M6T0_Exit_ts': M6T0_Exit_ts,          
          # Right hand enters left target
          'M6T1_Entry_ind' : M6T1_Entry_ind,
          'M6T1_Entry_ts' : M6T1_Entry_ts,
          # Right hand exits left target
          'M6T1_Exit_ind' : M6T1_Exit_ind,
          'M6T1_Exit_ts' : M6T1_Exit_ts,
          # Right hand enters center target
          'M6CT_Entry_ind' : M6CT_Entry_ind,
          'M6CT_Entry_ts' : M6CT_Entry_ts,
          # Right hand exits center target
          'M6CT_Exit_ind' : M6CT_Exit_ind,
          'M6CT_Exit_ts' : M6CT_Exit_ts,
          # Left hand enters right target
          'M7T0_Entry_ind' : M7T0_Entry_ind,
          'M7T0_Entry_ts' : M7T0_Entry_ts,
          # Left hand exits right target
          'M7T0_Exit_ind' : M7T0_Exit_ind,
          'M7T0_Exit_ts' : M7T0_Exit_ts,
          # Left hand enters left target
          'M7T1_Entry_ind' : M7T1_Entry_ind,
          'M7T1_Entry_ts' : M7T1_Entry_ts,
          # Left hand exits left target
          'M7T1_Exit_ind' : M7T1_Exit_ind,
          'M7T1_Exit_ts' : M7T1_Exit_ts,
          # Left hand enters center target
          'M7CT_Entry_ind' : M7CT_Entry_ind,
          'M7CT_Entry_ts' : M7CT_Entry_ts,          
          # Left hand exits center target
          'M7CT_Exit_ind' : M7CT_Exit_ind,
          'M7CT_Exit_ts' : M7CT_Exit_ts
         }
  
def EventComparator(tlist1, tlist2, tol_window):
  
  # This tests if each entry t_2 in tlist2 is within the range [t_1 + tol_lo, t_1 + tol_hi) 
  # for each entry in tlist1.  Output is an object array of size 
  # tlist1.size having each field contain a list of indices that point to all entries
  # in tlist2 that lie within the range. Note, for t + tol_lo to precede t, tol_lo must
  # be negative-valued.
  
  # Extract low tolerance and high tolerance boundaries from the tol_window tuple.
  (tol_lo, tol_hi) = tol_window
  
  # Generate index lists for subsequent iteration and filtering procedures.
  tlist2_ind = np.arange(0, tlist2.size)
  tlist1_ind = np.arange(0, tlist1.size)
  
  # Pre-allocate array dimensions for output. Note that the datatype is an
  # array of dicts.
  output = np.empty((tlist1.size,), dtype=dict)
  
  for i in tlist1_ind:
    
    filt = (tlist2 >= tlist1[i] + tol_lo) & (tlist2 < tlist1[i] + tol_hi)
    output[i]= {
                'within_tol_ts' : tlist2[filt], 
                'within_tol_ind' :tlist2_ind[filt]
               }
   
  return output

def RemoveEventsInTolWindowFiltGen(EventComparatorOutput):
    
    (nEvents,) = EventComparatorOutput.shape
    KeepEvent = (np.ones((nEvents,), dtype=int) == 1)
    
    for i in np.arange(0, nEvents):
        
        KeepEvent[EventComparatorOutput[i]['within_tol_ind']] = False
  
    return KeepEvent

def KeepOnlyEventsInTolWindowFiltGen(EventComparatorOutput):
    
    (nEvents,) = EventComparatorOutput.shape
    KeepEvent = (np.ones((nEvents,), dtype=int) == 0)
    
    for i in np.arange(0, nEvents):
        
        KeepEvent[EventComparatorOutput[i]['within_tol_ind']] = True
  
    return KeepEvent
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
import pandas as pd
import CalciumImagingBehaviorAnalysis as CIBA
from collections import defaultdict
from shutil import copyfile
from sys import exit
import sys

PathToAccessDir = '/home/thugwithyoyo/BoxFTP/MoxonLab/Projects/NHP Calcium Imaging/NHP_CalciumImaging/PlexonData/RetrackingData'
PathToTextFile = '/home/thugwithyoyo/Desktop/PathsToCSVs_201812.txt'
PathToWriteFile = '/home/thugwithyoyo/Desktop/EventCounts.xlsx'
PathToCacheDir = '/home/thugwithyoyo/Desktop/CacheDir'

ScriptDir = os.getcwd()

RefEventsList = [
                 'M6T0_Entry_ts', 'M6T0_Exit_ts', 
                 'M6T1_Entry_ts', 'M6T1_Exit_ts', 
                 'M6CT_Entry_ts', 'M6CT_Exit_ts',
                 'M7T0_Entry_ts', 'M7T0_Exit_ts',
                 'M7T1_Entry_ts', 'M7T1_Exit_ts',
                 'M7CT_Entry_ts', 'M7CT_Exit_ts'                
                ]

file = open(PathToTextFile, 'r')

FileList = file.readlines()

file.close()

# Count number of file entries in list.
nFiles = len(FileList)

# Remove junk characters from edges of strings.
for i in np.arange(0, nFiles):
    
    FileList[i] = FileList[i][2:-1]

#EventCountsByFileArray = np.empty((nFiles, 2*nChannelsToCount))
#EventCounts_df = pd.DataFrame(index=FileList, columns=RefEventsList)
EventCounts_Dict = defaultdict(dict)

for Filename in FileList:
    
    print('File to read: ' + Filename)
    # Construct path to file.
    FilePath = PathToAccessDir + '/' + Filename
    SavePath = PathToCacheDir + '/' + Filename.split(os.sep)[-1]
    
    # Download file to local cache
    # Navigate to cache directory
    os.chdir(PathToCacheDir)
    
    # adding exception handling
    try:
        
        copyfile(FilePath, SavePath)

    except IOError as e:

        print("Unable to copy file. %s" % e)
        exit(1)

    except:

        print("Unexpected error:", sys.exc_info())
        exit(1)
    
    print("\nFile copy done!\n")    
    
    #os.popen('cp ' + FilePath + ' .')
    #copyfile(FilePath, SavePath)
    
    # Run behavioral parser on current file.
    OutputDict = CIBA.CinePlexCSV_parser(SavePath)
    
    # Delete file from cache
    #os.popen('rm ' + Filename)
    try:
        
        os.remove(SavePath)
    
    except OSError as e:  ## if failed, report it back to the user ##
    
        print ("Error: %s - %s." % (e.filename, e.strerror))
    
    # Navigate back to current working directory
    os.chdir(ScriptDir)
    
    # Count number of timestamps in each element of OutputDict.
    for el in RefEventsList:

        # Write counts to corresponding columns in dataframe. The row 
        # will be given the name of the file. 
        #EventCounts_df.set_value(Filename, el, OutputDict[el].size)
        #EventCounts_df.at[Filename, el] = OutputDict[el].size
        EventCounts_Dict[Filename][el] = OutputDict[el].size
        
EventCounts_df = pd.DataFrame(data=EventCounts_Dict)
EventCounts_df = (EventCounts_df.T)
EventCounts_df.to_excel(PathToWriteFile, columns=RefEventsList.sort())
    
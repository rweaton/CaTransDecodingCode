#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 08:09:15 2019

@author: thugwithyoyo
"""
import os
import shelve

# This script requires that the path variable SavePath be pre-defined
# in the scope from which this script is called.
#
# Example: code calling this script should be as follows:
#
#    SavePath = '/home/thugwithyoyo/CaTransDecoding/Output/TestSave2'
#    exec(open('./ShelveWorkspaceScript.py').read())

# Save path of directory that contains this script
ScriptDir = os.getcwd()
    
# Use path from function argument to change to save directory and write file
# with name contained in Filename.
(PathToFile, Filename) = os.path.split(SavePath)

# Change to directory that will contain file to be saved.
os.chdir(PathToFile)
 
# Open a shelf object to contain workspace variables.
my_shelf = shelve.open(Filename)    
#PathToSaveFile = root.filename
#my_shelf = shelve.open(PathToSaveFile, 'n') # 'n' for new

# Iterate through list of "global" variables and write to file. Write each
# element to the my_shelf dictionary object. 
#
# Remember that globals() actually generates a list of the variables restricted 
# to within scope of the calling function.
for key in dir():
    
    try:
        
        #my_shelf[key] = globals()[key]
        my_shelf[key] = locals()[key]
        
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))

# Close shelf object after variables have been written to file
my_shelf.close()

# Return to directory that contains this script.
os.chdir(ScriptDir)
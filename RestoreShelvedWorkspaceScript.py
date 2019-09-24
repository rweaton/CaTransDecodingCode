#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 07:46:36 2019

@author: thugwithyoyo
"""
import os
import shelve

# This script requires that the path variable RestoreFilePath be pre-defined
# in the scope from which this script is called.
#
# Example: code calling this script should be as follows:
#
#    RestoreFilePath = '/home/thugwithyoyo/CaTransDecoding/Output/TestSave2.dat'
#    exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Save path of directory that contains this script
ScriptDir = os.getcwd()

(PathToFile, Filename) = os.path.split(RestoreFilePath)

# Change to directory that contains file to be loaded.
os.chdir(PathToFile)
 
# Open a shelf dictionary object that contains stored variables
my_shelf = shelve.open(os.path.splitext(Filename)[0])

# Iterate through dictionary contents and write each to the global variables
# list.
for key in my_shelf:
    
    globals()[key]=my_shelf[key]

# Close the shelve object
my_shelf.close()

# Return to script directory
os.chdir(ScriptDir)
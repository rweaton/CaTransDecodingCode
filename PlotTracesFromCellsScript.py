#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 09:03:29 2019

@author: thugwithyoyo
"""

# Be sure to run PlotPeriEventAverages_TuningSorted (with resultant workspace
# loaded) prior to running this script.

# 2018-12-18-11-20-21 cells to select
#CellsToSelectList = ['C031', 'C085', 'C025', 'C081', 'C099', 'C101']
CellsToSelectList=['C025', 'C031', 'C085']

# select cells
CellsToSelectList = np.array(CellsToSelectList)
(NumCellsToSelect,) = CellsToSelectList.shape


# Generate averaged z-scored traces version of plot
fig, axes  = plt.subplots(1, NumCellsToSelect)

fig.suptitle('Averaged Z-SCORED $\Delta$F trace averages from cells')

for i in np.arange(0, NumCellsToSelect):
    
    axes[i].set_title(str(CellsToSelectList[i]))
    
    if i == 0:
        axes[i].set_ylabel('Average z-scored $\Delta$F')
    
        axes[i].set_xlabel('time relative to zone-entry (sec.)')
    
    for j in np.arange(0, NumEvents):
        
        axes[i].plot(RelTimeVec, 
            AveragedTracesMatrices[SortProcessingDict['SortIndices']\
                                   [Tuning_Frame['HeatMapRowIndex'].loc[CellsToSelectList[i]]], :, j])

# Generate non-z-scored, raw trace averages version
    
# Reload orignal cell fluorescence data frame
#RestoreFilePath = SavePath +'.dat'
root = tk.Tk()
RestoreFilePath = askopenfilename()
root.withdraw()

# Determine parent directory and filename from complete path.
drive, path_and_file = os.path.splitdrive(RestoreFilePath)
path, file = os.path.split(path_and_file)

exec(open('./RestoreShelvedWorkspaceScript.py').read())

# Re-initialize settings
### Specify a dictionary of parameters for tuning calcluation and plotting.
# Time domain to consider for calculation of tuning index
ParamsDict['SearchDomain'] = [0., 1.]

# Specify the method to be used to calculate tuning method. 
#ParamsDict['TuningType'] = 'WeightedDifference'
ParamsDict['TuningType'] = 'DiffOfAvgsOverSumOfMagOfAvgs'

# For the WeightedDifference method, specify the standard deviation of the 
# Gaussian weighting function.
ParamsDict['StdDev'] = 0.75

# Specify the number of bins to use to generate the histogram.
ParamsDict['NumHistoBins'] = 10

# Specify the list of quantiles to plot as vertical lines on histogram
ParamsDict['QuantilesList'] = [0.25, 0.75]

# The magnitude of the tuning index value that separates tuned cells from 
# non-tuned cells.
ParamsDict['TuningCutoffLevel'] = 0.50

# Name of the color scheme used for plotting heat map.
ColorMap = 'seismic'

RefEventsList = RefEventsDict['RefEventsList']
SortingRefEvent = RefEventsList[1]

nEvents = len(RefEventsList)

# Acquire number of columns in fluorescence trace dataframe
nCols = CellFluorTraces_Frame.shape[1]

# BYPASS zScoring the calcium traces prior to peri-event extraction
# Z-scoring and replacing the z-scored traces in the array takes place here.

# Recompile PeriEvent activity dict to include Full Domain in boundary window.
PeriEventExtractorDict = PeriEventExtractor_Trace(BehavDict, 
                                                  CellFluorTraces_Frame, 
                                                  RefEventsDict, 
                                                  ParamsDict['BoundaryWindow'])

# Count the number of trials and determine length of each row in the formatted
# Peri-event activity array.
NumTrials, NumCellsByNumSamples = PeriEventExtractorDict['PEA_Array'].shape

# Count the total number of columns in the fluorescence trace dataframe
(_, NumColumns) = CellFluorTraces_Frame.shape

# 
NumCells = NumColumns - 1

NumSamples = int(NumCellsByNumSamples/NumCells)

RelTimeVec = np.linspace(ParamsDict['BoundaryWindow'][0], 
                         ParamsDict['BoundaryWindow'][1], 
                         num=NumSamples, endpoint=False)

NumEvents = len(RefEventsList)

AveragedTracesMatrices = np.empty((NumCells, NumSamples, NumEvents))

i = 0
for RefEvent in RefEventsList:
    
    TracesFilt = PeriEventExtractorDict['TrialIndicesByEventDict'][RefEvent]
    
    AveragedTraces = np.mean(PeriEventExtractorDict['PEA_Array'][TracesFilt,:],
                                      axis=0)
    
    
    AveragedTracesMatrices[:,:,i] = np.reshape(AveragedTraces, (NumCells, NumSamples))
    
    i += 1

fig, axes  = plt.subplots(1, NumCellsToSelect)

fig.suptitle('Averaged (raw) $\Delta$F trace averages from cells')

for i in np.arange(0, NumCellsToSelect):
    
    axes[i].set_title(str(CellsToSelectList[i]))
    
    if i == 0:
        axes[i].set_ylabel('Average $\Delta$F')
    
        axes[i].set_xlabel('time relative to zone-entry (sec.)')
        
    for j in np.arange(0, NumEvents):
        
        axes[i].plot(RelTimeVec, 
            AveragedTracesMatrices[SortProcessingDict['SortIndices']\
                                   [Tuning_Frame['HeatMapRowIndex'].loc[CellsToSelectList[i]]], :, j])
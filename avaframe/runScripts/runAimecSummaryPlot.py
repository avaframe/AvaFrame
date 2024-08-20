"""
    Run script for creating aimec summary plot
    distribution of scalar measures
"""
# Load modules
# importing general python modules
import pathlib
import pandas as pd
import numpy as np
import logging

# Local imports
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.out3Plot import outAIMEC
from avaframe.in3Utils import logUtils

# ++++++++++user input+++++++
# fetch input resAnalysis files from AIMEC analysis
# provide directory where aimec *_resAnalysisDF.csv files are located
resAnalysisDir = ''
# set paths to aimec resAnalysisDF csv files to be loaded and then compared in the summary plots
df1File = pathlib.Path(resAnalysisDir, 'FILENAME1')
df2File = pathlib.Path(resAnalysisDir, 'FILENAME2')
# set a directory where new summary plots shall be saved
resultPath = pathlib.Path(resAnalysisDir, 'plots')
# name for test cases that describe the individual aimec analysis files (e.g. two different avalanche directoy names)
testCase1 = 'avaTest1'
testCase2 = 'avaTest2'
# +++++++++++++++++++++++++++

# result variables that are analysed
resultVariables = ['sRunout', 'runoutAngle', 'pfvFieldMax', 'pftFieldMax']

# create local logger
log = logging.getLogger(__name__)
logName = 'runAimecSummaryPlot'

# Start logging
log = logUtils.initiateLogger(resAnalysisDir, logName)
log.info('MAIN SCRIPT')
log.info('Results saved to %s', resultPath)

# read dataframes from all test cases
df1 = pd.read_csv(df1File)
df2 = pd.read_csv(df2File)

# add testCase name and number
df1['testCase'] = [testCase1] * len(df1)
df2['testCase'] = [testCase2] * len(df2)

# concatenate dataframes from all test cases
dfAll = pd.concat([df1, df2], ignore_index=True)
dfPoints = pd.DataFrame(data={'testCasePoints': [testCase1+'_points'] * len(df1) + [testCase2+'_points'] * len(df2)})
dfAll = pd.concat([dfAll, dfPoints], axis=1)

# compute interquartile range
statisticsDict = {}
for testC in [testCase1, testCase2]:
    for resVar in resultVariables:
        dfAnalyze = dfAll[dfAll['testCase'] == testC]
        resultField = dfAnalyze[resVar].to_numpy()
        resultField = resultField[np.isnan(resultField) == False]
        q1, q3 = np.percentile(resultField, [25, 75])
        median = np.median(resultField)
        log.info('testCase: %s, result variable: %s' % (testC, resVar))
        log.info('q1: %.2f, q3: %.2f, iqr: %.2f, median: %.2f' % (q1, q3, q3-q1, np.median(resultField)))
        iqrpercent = (100. / median) * (q3-q1)
        log.info('iqr percent of median %.02f' % iqrpercent)
        statisticsDict[testC] = {resVar: {'q3': q3, 'q1': q1, 'iqr': q3-q1, 'median': median}}


# create figure and save
pathDict = {'pathResult': resultPath}
fU.makeADir(resultPath)
orderList = [testCase1, testCase1+'_points', testCase2, testCase2+'_points']
outAIMEC.boxScalarMeasures(pathDict, dfAll, name='box', orderList=orderList)
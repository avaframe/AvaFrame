"""
    Pytest for module com1DFA
"""

#  Load modules
import configparser
import pathlib
import shutil


from avaframe.com5SnowSlide import com5SnowSlide
from avaframe.in3Utils import cfgUtils


def test_snowSlide(tmp_path):
    """ test calling of snowSlide """

    testDir = pathlib.Path(__file__).parents[0]
    inputDir = testDir / 'data' / 'testCom1DFA'
    avaDir = pathlib.Path(tmp_path, 'testCom1DFA')
    shutil.copytree(inputDir, avaDir)

    # setup requuired input data
    cfgMain = configparser.ConfigParser()
    cfgMain['FLAGS'] = {'showPlot': 'False', 'savePlot': 'True', 'ReportDir': 'True', 'reportOneFile': 'True',
        'debugPlot': 'False'}
    cfgMain['MAIN'] = {'avalancheDir': str(avaDir), 'nCPU': '1'}

    snowSlideCfg = cfgUtils.getModuleConfig(com5SnowSlide)
    snowSlideCfg['com1DFA_com1DFA_override']['dt'] = '1'
    snowSlideCfg['com1DFA_com1DFA_override']['meshCellSize'] = '5'
    snowSlideCfg['com1DFA_com1DFA_override']['tEnd'] = '5'
    snowSlideCfg['com1DFA_com1DFA_override']['entThFromShp'] = 'False'
    snowSlideCfg['com1DFA_com1DFA_override']['entTh'] = '0.3'

    # call function to be tested
    _, _, _, simDF = com5SnowSlide.com5SnowSlideMain(cfgMain, snowSlideCfg)

    # print('simDF', simDF.to_string())
    assert simDF['tEnd'].iloc[0] == 5
    assert simDF['relThFromShp'].iloc[0] == 'True'
    assert simDF['snowSlide'].iloc[0] == 1

"""Test functions for in4Region module."""

import pathlib
import configparser
import shutil
from avaframe.in4Region import splitInputs

def test_splitInputsMain(tmp_path):
    """Test splitInputsMain function using pre-generated test data."""
    # Set up test data
    test_data_dir = pathlib.Path(__file__).parent / 'data' / 'testIn4Region'
    inputDir = tmp_path / 'input'
    shutil.copytree(test_data_dir, inputDir)
    outputDir = tmp_path / 'output'
    
    # Configure test parameters
    cfg = configparser.ConfigParser()
    cfg['GENERAL'] = {'bufferSize': '50'}
    cfgMain = configparser.ConfigParser()
    cfgMain['FLAGS'] = {'createReport': 'True', 'savePlot': 'True'}
    
    # Run function
    splitInputs.splitInputsMain(inputDir, outputDir, cfg, cfgMain)
    
    # Verify outputs
    assert outputDir.exists()
    
    # Check group directories
    groupDirs = list(outputDir.glob('group*'))
    assert len(groupDirs) == 2
    
    for groupDir in groupDirs:
        # Check directory structure
        assert (groupDir / 'Inputs').exists()
        assert (groupDir / 'Inputs' / 'REL').exists()
        assert len(list((groupDir / 'Inputs').glob('*.tif'))) == 1
        
        # Check release areas were split by scenarios
        relDir = groupDir / 'Inputs' / 'REL'
        assert len(list(relDir.glob('*.shp'))) == 2
        
        # Check optional inputs
        assert (groupDir / 'Inputs' / 'ENT').exists()
        assert len(list((groupDir / 'Inputs' / 'ENT').glob('*.shp'))) == 1
        assert (groupDir / 'Inputs' / 'RES').exists()
        assert len(list((groupDir / 'Inputs' / 'RES').glob('*.shp'))) == 1
    
    # Check reports
    assert (outputDir / 'splitInputs_scenarioReport.txt').exists()
    assert (outputDir / 'splitInputs_visualReport_basic.png').exists()
    assert (outputDir / 'splitInputs_visualReport_optional.png').exists()

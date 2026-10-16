#!/usr/bin/env encoding
# -*- coding: utf-8 -*-

"""
    Run script to execute com4FlowPy benchmarks using local custom datasets
    defined in local config files (avaframeCfg.ini and com4FlowPyCfg.ini).
    Profiles realTime (wall-time) and peak RAM usage (in GiB) utilizing all
    available CPU cores. Reports spatial DEM properties and active global inputs.
    Produces a Markdown report under tests/reportsCom4FlowPy/
"""

import time
import pathlib
import numpy as np
from datetime import datetime
import os
import threading
import psutil
import logging
import shutil

# Local imports
import avaframe as avaf
from avaframe.com4FlowPy import com4FlowPy
from avaframe.runCom4FlowPy import readFlowPyinputs
from avaframe.in3Utils import fileHandlerUtils as fU
from avaframe.in3Utils import initializeProject as initProj
from avaframe.in3Utils import cfgUtils
from avaframe.in3Utils import logUtils
import avaframe.in2Trans.rasterUtils as rasterUtils


class PeakMemoryTrackerGiB:
    """Recursively tracks memory of parent and all child processes in GiB."""
    def __init__(self, interval=0.01):
        self.interval = interval
        self.peakMemGiB = 0.0
        self.running = False
        self._thread = None

    def _track(self):
        parent = psutil.Process()
        while self.running:
            try:
                # Add memory of parent process
                totalMem = parent.memory_info().rss
                # Add memory of all children recursively
                for child in parent.children(recursive=True):
                    try:
                        totalMem += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                totalMemGiB = totalMem / (1024 * 1024 * 1024)
                if totalMemGiB > self.peakMemGiB:
                    self.peakMemGiB = totalMemGiB
            except Exception:
                pass
            time.sleep(self.interval)

    def start(self):
        self.peakMemGiB = 0.0
        self.running = True
        self._thread = threading.Thread(target=self._track, daemon=True)
        self._thread.start()

    def stop(self):
        self.running = False
        if self._thread:
            self._thread.join(timeout=1.0)
        return self.peakMemGiB


def main():
    # avaframe directory
    _avaframeDir = pathlib.Path(avaf.__file__).parents[0]
    outDir = _avaframeDir / 'tests' / 'reportsCom4FlowPy'
    fU.makeADir(outDir)

    logName = 'runCustomBenchmarkRAMcom4FlowPy'
    log = logUtils.initiateLogger(outDir, logName)
    log.info('Starting Custom Com4FlowPy RAM & Performance Benchmark')

    # Load local configurations
    cfgMain = cfgUtils.getGeneralConfig()
    cfg = cfgUtils.getModuleConfig(com4FlowPy)
    cfgGen = cfg["GENERAL"]

    useCustomPaths = cfg.getboolean("PATHS", "useCustomPaths", fallback=False)

    if useCustomPaths:
        log.info("useCustomPaths is True - utilizing custom paths defined in com4FlowPyCfg.ini [PATHS]")
        cfgCustomPaths = cfg["PATHS"]
        
        cfgPath = {}
        workDir = pathlib.Path(cfgCustomPaths["workDir"])
        res_dir = workDir / f"res_{cfgUtils.cfgHash(cfg)}"
        temp_dir = res_dir / "temp"
        
        cfgPath["workDir"] = workDir
        cfgPath["outDir"] = res_dir
        cfgPath["resDir"] = res_dir
        cfgPath["tempDir"] = temp_dir
        cfgPath["demPath"] = pathlib.Path(cfgCustomPaths["demPath"])
        cfgPath["releasePath"] = pathlib.Path(cfgCustomPaths["releasePath"])
        cfgPath["relIdPath"] = pathlib.Path(cfgCustomPaths["relIdPath"]) if cfgCustomPaths.get("relIdPath") else pathlib.Path("")
        cfgPath["infraPath"] = pathlib.Path(cfgCustomPaths["infraPath"]) if cfgCustomPaths.get("infraPath") else pathlib.Path("")
        cfgPath["forestPath"] = pathlib.Path(cfgCustomPaths["forestPath"]) if cfgCustomPaths.get("forestPath") else pathlib.Path("")
        cfgPath["varUmaxPath"] = pathlib.Path(cfgCustomPaths["varUmaxPath"]) if cfgCustomPaths.get("varUmaxPath") else pathlib.Path("")
        cfgPath["varAlphaPath"] = pathlib.Path(cfgCustomPaths["varAlphaPath"]) if cfgCustomPaths.get("varAlphaPath") else pathlib.Path("")
        cfgPath["varExponentPath"] = pathlib.Path(cfgCustomPaths["varExponentPath"]) if cfgCustomPaths.get("varExponentPath") else pathlib.Path("")
        cfgPath["deleteTemp"] = cfgCustomPaths.get("deleteTempFolder", "False")
        cfgPath["outputFileFormat"] = cfgCustomPaths.get("outputFileFormat", ".tif")
        cfgPath["outputFiles"] = cfgCustomPaths.get("outputFiles", "zDelta")
        cfgPath["outputNoDataValue"] = cfgCustomPaths.getfloat("outputNoDataValue", fallback=-9999.0)
        cfgPath["useCompression"] = cfgCustomPaths.getboolean("useCompression", fallback=False)
        cfgPath["customDirs"] = "True"
        
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
    else:
        log.info("useCustomPaths is False - utilizing AvaFrame directory structure from avaframeCfg.ini")
        # Retrieve paths to custom data
        avalancheDir = cfgMain["MAIN"]["avalancheDir"]
        cfgPath = readFlowPyinputs(avalancheDir, cfg, log)
        # Clean input directory of old work and output files
        initProj.cleanSingleAvaDir(avalancheDir, deleteOutput=False)
        cfgPath["customDirs"] = "False"

    demPath = pathlib.Path(cfgPath["demPath"])
    releasePath = pathlib.Path(cfgPath["releasePath"])

    # 1. Analyze DEM Spatial Properties & Size on Disk
    log.info(f"Analyzing DEM: {demPath}")
    header = rasterUtils.readRasterHeader(demPath)
    
    ncols = header["ncols"]
    nrows = header["nrows"]
    cellsize = header["cellsize"]
    total_cells = ncols * nrows
    
    x_extent = ncols * cellsize
    y_extent = nrows * cellsize
    
    dem_disk_bytes = os.path.getsize(demPath)
    dem_disk_mb = dem_disk_bytes / (1024 * 1024)

    # 2. Check Global/Permanent Inputs Status
    infra_active = cfgGen.getboolean("infra")
    forest_active = cfgGen.getboolean("forest")
    var_umax_active = cfgGen.getboolean("variableUmaxLim")
    var_alpha_active = cfgGen.getboolean("variableAlpha")
    var_exponent_active = cfgGen.getboolean("variableExponent")

    # 3. Setup CPU Cores for Multiprocessing
    maxCores = psutil.cpu_count(logical=True) or 1
    cfgGen["cpuCount"] = str(maxCores)
    cfgGen["overwriteResults"] = "reRunAndOverwrite"

    # Set outputs matching custom configs
    cfgPath["tempDir"] = cfgPath["workDir"] / "temp" if not useCustomPaths else cfgPath["tempDir"]
    fU.makeADir(cfgPath["tempDir"])
    cfgPath["deleteTemp"] = "False"
    cfgPath["outputFiles"] = cfg["PATHS"]["outputFiles"]
    cfgPath["outputNoDataValue"] = cfg["PATHS"].getfloat("outputNoDataValue")
    cfgPath["useCompression"] = cfg["PATHS"].getboolean("useCompression")
    cfgPath["uid"] = cfgUtils.cfgHash(cfg)
    cfgPath["timeString"] = datetime.now().strftime("%Y%m%d_%H%M%S")

    # 4. Profile Memory and Execution Time
    log.info("Starting Custom Simulation and Tracking Memory...")
    tracker = PeakMemoryTrackerGiB(interval=0.01)
    tracker.start()

    _startDate = datetime.now()
    startTime = time.perf_counter()
    com4FlowPy.com4FlowPyMain(cfgPath, cfgGen)
    endTime = time.perf_counter()
    _endDate = datetime.now()

    peakMemGiB = tracker.stop()
    wallTime = endTime - startTime

    log.info(f"Custom Benchmark finished. Wall-Time: {wallTime:.2f}s, Peak RAM: {peakMemGiB:.3f} GiB")

    # 5. Write Report to tests/reportsCom4FlowPy/customBenchmarkRAMcom4FlowPyReport.md
    reportFile = outDir / 'customBenchmarkRAMcom4FlowPyReport.md'
    
    with open(reportFile, 'w') as pfile:
        pfile.write('# Custom Com4FlowPy Performance and RAM Benchmark Report \n\n')
        pfile.write('Profiles the peak RAM (in GiB) and wall-time of com4FlowPy using custom dataset and local INI configurations.\n\n')
        pfile.write('* * * \n')
        pfile.write(f'**tested AvaFrame Version:** {avaf.version.getVersion()}  \n')
        pfile.write(f'**CPU Cores Utilized:** {maxCores} &check;  \n')
        pfile.write(f'**Run Started:** {_startDate.strftime("%Y-%m-%d %H:%M:%S")}  \n')
        pfile.write(f'**Run Finished:** {_endDate.strftime("%Y-%m-%d %H:%M:%S")}  \n')
        pfile.write('* * * \n\n')

        pfile.write('## 🚀 Custom Run Performance Metrics\n\n')
        pfile.write(f'* **Wall-Time (Elapsed Duration):** `{wallTime:.3f} seconds`  \n')
        pfile.write(f'* **Peak RAM Usage:** `{peakMemGiB:.3f} GiB`  \n')
        pfile.write(f'* **Simulation UID:** `{cfgPath["uid"]}`  \n\n')

        pfile.write('## 🗺️ Spatial DEM Metadata\n\n')
        pfile.write(f'* **DEM Path:** `{demPath}`  \n')
        pfile.write(f'* **Spatial Dimensions (N Cols x N Rows):** `{ncols} x {nrows}` cells  \n')
        pfile.write(f'* **Total Cells ($N$):** `{total_cells:,}` cells  \n')
        pfile.write(f'* **Resolution (Cell Size):** `{cellsize} meters`  \n')
        pfile.write(f'* **Horizontal Extent (x-axis):** `{x_extent:,.2f} meters`  \n')
        pfile.write(f'* **Vertical Extent (y-axis):** `{y_extent:,.2f} meters`  \n')
        pfile.write(f'* **File Size on Disk:** `{dem_disk_mb:.2f} MB`  \n\n')

        pfile.write('## ⚙️ Active Global Inputs & Configs\n\n')
        pfile.write('| Config Input | Status | Details / Path |\n')
        pfile.write('| :--- | :---: | :--- |\n')
        pfile.write(f"| **Infrastructure** | {'Enabled &check;' if infra_active else 'Disabled &cross;'} | `{cfgPath['infraPath'] if infra_active else ''}` |\n")
        pfile.write(f"| **Forest resistance** | {'Enabled &check;' if forest_active else 'Disabled &cross;'} | `{cfgPath['forestPath'] if forest_active else ''}` |\n")
        pfile.write(f"| **Variable Umax Limit** | {'Enabled &check;' if var_umax_active else 'Disabled &cross;'} | `{cfgPath['varUmaxPath'] if var_umax_active else ''}` |\n")
        pfile.write(f"| **Variable Alpha** | {'Enabled &check;' if var_alpha_active else 'Disabled &cross;'} | `{cfgPath['varAlphaPath'] if var_alpha_active else ''}` |\n")
        pfile.write(f"| **Variable Exponent** | {'Enabled &check;' if var_exponent_active else 'Disabled &cross;'} | `{cfgPath['varExponentPath'] if var_exponent_active else ''}` |\n")
        pfile.write('\n* * * \n')

    log.info(f"Custom benchmarking complete. Report written to {reportFile}")


if __name__ == "__main__":
    main()

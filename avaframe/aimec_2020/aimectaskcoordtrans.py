#! /usr/bin/python

from aimectask import *

import re
import aimecdata
import aimecresult
import esrigrid
import numpy as np
import math
import time
import os


from coordtrans_process import processData
from coordtrans_process import analyze
from coordtrans_process import aimecOut


# from PyQt5 import QtCore


class AimecTaskCoordTrans(AimecTask):

    def __init__self(self):

        super(AimecTaskCoordTrans, self).__init__()

    def cliName(self):

        return '-trans'

    def name(self):

        return 'Coordinate Transformation'

    def resultPath(self):

        return 'coordtrans'

    def description(self):

        return ('Simualtion result transformation from global coordinate system to a path dependend one.')

    def validateData(self, data):

        # TODO fixme
        return True

    def run(self, data, callBack=None):

        startTime = time.time()

        # -----------------------------------------------------------
        # create new raster + preparing new raster assignment function
        # -----------------------------------------------------------

        pressurefileList = [str(data.pathPressure) +
                            '/' +
                            str(name) for name in
                            sorted(os.listdir(data.pathPressure)) if os.path.isfile(os.path.join(data.pathPressure, name))]

        set_name = pressurefileList[0].split('/')[-3]
        project_name = str(data.pathAvalanchePath).split('/')[-4]
        path_name = str(data.pathAvalanchePath).split('/')[-1]

        deskewedRasterInd = processData.processDataInd(pressurefileList,
                                                       str(data.pathAvalanchePath),
                                                       data.domainWidth,
                                                       data.with__(data.pathDepoArea),
                                                       str(data.pathDepoArea),
                                                       data.with__(data.pathDocDamage),
                                                       str(data.pathDocDamage),
                                                       True,  # str(data.pathMapResult),
                                                       str(data.pathResult),
                                                       data.with__(data.pathResult),
                                                       data.calcPressureLimit)

        # -----------------------------------------------------------
        # transform pressure_data and depth_data in new raster
        # -----------------------------------------------------------

        deskewedRasterPressure = processData.assignData(pressurefileList, deskewedRasterInd)

        if data.with__(data.pathFlowHeight):
            print('[MAIN] ----------------------------------')
            print('[MAIN] assigning depth data in new raster')
            # assign depth data

            depthfileList = [str(data.pathFlowHeight) +
                                '/' +
                                str(name) for name in
                                sorted(os.listdir(data.pathFlowHeight)) if os.path.isfile(os.path.join(data.pathFlowHeight, name))]

            deskewedRasterDepth = processData.assignData(depthfileList, deskewedRasterInd)
        else:
            deskewedRasterDepth = None

        if data.with__(data.pathDHM):
            # if dhm delta h analysis
            # Achtung Fehler in SamosAT: Druckraster und DHM-Raster stimmen nicht exakt ueberein!
            # Eventuell shift in assignData beruecksichtigen
            deskewedRasterDHM = processData.assignData([str(data.pathDHM)], deskewedRasterInd)
            dhm_name = str(data.pathDHM).split('/')[-1]

        else:
            deskewedRasterDHM = None
            dhm_name = None

        # -----------------------------------------------------------
        # create doku_mask
        # -----------------------------------------------------------

        if data.with__(data.pathDepoArea):
            print('[MAIN] ----------------------------------')
            print('[MAIN] create and transferring doku in new raster')
            deskewedRasterDoku = processData.createMask(deskewedRasterInd,
                                                        str(data.pathDepoArea),
                                                        data.calcPressureLimit)
            doku_name = str(data.pathDepoArea).split('/')[-1]
        else:
            deskewedRasterDoku = None
            doku_name = str(pressurefileList[0]).split('/')[-1]

        # -----------------------------------------------------------
        # create damages_mask
        # -----------------------------------------------------------
        if data.with__(data.pathDocDamage):
            print('[MAIN] ----------------------------------')
            print('[MAIN] create and transferring damages in new raster')
            deskewedRasterDamage = processData.createMask(deskewedRasterInd,
                                                          str(data.pathDocDamage),
                                                          data.calcPressureLimit)
            damage_name = str(data.pathDocDamage).split('/')[-1]
        else:
            deskewedRasterDamage = None
            damage_name = ''

        # -----------------------------------------------------------
        # analyze
        # -----------------------------------------------------------

        # analyze doku
        print('[MAIN] analyzing dokumentation data')
        doku, runout_doku, delta_h, elevRel = analyze.analyzeDocu(data.calcPressureLimit,
                                                                  pressurefileList,
                                                                  deskewedRasterInd,
                                                                  deskewedRasterPressure,
                                                                  deskewedRasterDepth,
                                                                  deskewedRasterDoku,
                                                                  deskewedRasterDHM,
                                                                  data.with__(data.pathFlowHeight),
                                                                  data.with__(data.pathDHM),
                                                                  data.with__(data.pathDepoArea)
                                                                  )

        # analyze damages
        if data.with__(data.pathDocDamage):
            dam_mean, dam_max = analyze.analyzeImpact(pressurefileList, deskewedRasterPressure,
                                                      deskewedRasterDepth,
                                                      deskewedRasterDamage,
                                                      data.with__(data.pathFlowHeight))
        else:
            dam_mean = np.zeros((1, len(pressurefileList)))
            dam_max = np.zeros((1, len(pressurefileList)))

        # analyze mass / entrainment
        if data.with__(data.pathMass):
            print('[MAIN] analyzing entrainment data')
            # determine growth index from entrainment data
            [relMass,
             entMass,
             gr_index,
             gr_grad] = analyze.analyzeEntrainmentdata(pressurefileList)
            # gr_index [initial mass, max mass / initial mass]
        else:
            relMass = 0.
            entMass = np.zeros((len(pressurefileList)))
            gr_index = np.zeros((len(pressurefileList)))
            gr_grad = np.zeros((len(pressurefileList)))

        # analyze pressure_data and depth_data
        # determine runount, AMPP, AMD, FS,
        print('[MAIN] analyzing data in path coordinate system')
        [runout,
         runout_mean,
         AMPP,
         MMPP,
         AMD,
         MMD] = analyze.analyzeDataWithDepth(deskewedRasterInd,
                                             data.calcPressureLimit,
                                             pressurefileList, deskewedRasterPressure,
                                             deskewedRasterDepth,
                                             False,  # data.with__(data.pathMapResult),
                                             str(data.pathResult),
                                             data.with__(data.pathResult))

        # -----------------------------------------------------------
        # result visualisation + report
        # -----------------------------------------------------------
        aimecOut.result_visu(pressurefileList, str(data.pathAvalanchePath),
                             runout, AMPP, doku, gr_index, data.calcPressureLimit,
                             str(data.pathResult))

        # -----------------------------------------------------------
        # write results to file
        # -----------------------------------------------------------
        if data.with__(data.pathResult):
            #     write output file for postprocessing
            out_header = ''.join(['project_name: ',  str(project_name), '\n',
                                  'set_name: ', str(set_name), '\n',
                                  'path: ', str(path_name), '\n',
                                  'docu: ', str(doku_name), '\n',
                                  'dhm: ', str(dhm_name), '\n',
                                  'domain_width: ', str(data.domainWidth), '\n',
                                  'pressure_limit: ', str(data.calcPressureLimit), '\n',
                                  'runout_doku: ', str(runout_doku), '\n',
                                  'fall_height: ', str(delta_h), '\n',
                                  'release_mass: ', str(relMass), '\n',
                                  'elevation_release: ', str(elevRel),  '\n',
                                  'filenr, runout, AMPP, MMPP, entMass, growth_index, AMD, MMD, TPs, FNs, FPs, TNs, TP_depth, TP_pressure, damages_mean (%i)\n' % len(dam_mean)])
            outname = ''.join([str(data.pathResult), '/', self.resultPath(), '/',
                               str(set_name), '_pl', str(int(data.calcPressureLimit)), '_w', str(int(data.domainWidth)), '.txt'])

            print('[MAIN] write output file: %s' % outname)
            resfile = [runout, AMPP, MMPP, entMass, gr_index, AMD, MMD,
                       doku[0], doku[1], doku[2], doku[3], doku[4], doku[6]]
            resfile.extend(dam_mean)
            aimecOut.result_write(pressurefileList, resfile, outname, out_header)
            print('[MAIN] --- all done ---')

        # -----------------------------------------------------------
        # analyze numeric data and write to file
        # -----------------------------------------------------------
        if data.with__(data.pathNumInfo):
            numvar_results = np.zeros((4, len(pressurefileList)))
            numfileList = [str(data.pathNumInfo) +
                                '/' +
                                str(name) for name in
                                sorted(os.listdir(data.pathNumInfo)) if os.path.isfile(os.path.join(data.pathNumInfo, name))]

            for i in range(len(numfileList)):
                # load data
                numvar_results[:, i] = np.loadtxt(numfileList[i], skiprows=1)
                # total_mass, particles, mpp, total_time
            if data.with__(data.pathResult):
                out_header_num = ''.join(['project_name: ',  str(project_name), '\n',
                                          'parameter_set: ', str(set_name), '\n',
                                          'path: ', str(path_name), '\n',
                                          'docu: ', str(doku_name), '\n',
                                          'dhm: ', str(dhm_name), '\n',
                                          'domain_width: ', str(data.domainWidth), '\n',
                                          'pressure_limit: ', str(data.calcPressureLimit), '\n',
                                          'runout_doku: ', str(runout_doku), '\n',
                                          'fall_height: ', str(delta_h), '\n',
                                          'release_mass: ', str(relMass), '\n',
                                          'elevation_release: ', str(elevRel),  '\n',
                                          'filenr, total_mass, particles, mpp, total_time\n'])
                outname_num = ''.join([str(data.pathResult), '/', self.resultPath(), '/',
                                       str(set_name), '_numresults.txt'])
                print('[MAIN] write output file: %s' % outname_num)
                aimecOut.result_write(pressurefileList, numvar_results, outname_num, out_header_num)
                print('[MAIN] --- all done ---')

        # -----------------------------------------------------------
        # FIN --> post-processing
        # -----------------------------------------------------------
        print('[MAIN] --- finished --- ')

        endTime = time.time()

        print(('Took %s seconds to calculate.' % (endTime - startTime)))

        result = aimecresult.AimecResult(self, aimecresult.AimecReference([runout_doku, float('nan'), float('nan'), float('nan'), float('nan')],
                                                                          ['runout', 'tp', 'fn', 'fp', 'tn'], [1, 2, 3, 4, 5]))
        i = 0
        for r in runout:
            result.simulations.append(aimecresult.AimecResultSingle(i, [runout[i], doku[0][i], doku[1][i], doku[2][i], doku[3][i]],
                                                                    [math.fabs(runout_doku-runout[i]), doku[0][i], doku[1][i], doku[2][i], doku[3][i]]))
            i += 1

        result.save(data.pathResult, self.resultPath())

        return result


if __name__ == '__main__':

    infile = '/home/matthiastonnel/Documents/github/Avaframe/avaframe/aimec_2020/NameOfAvalanche/Outputs/aimecPathFile.txt'
    f = open(infile, 'r')
    indata = aimecdata.fromString(f.read())

    task = AimecTaskCoordTrans()

    print('data is ok = ' + task.validateData(indata).__str__())

    print(task.run(indata).__str__())

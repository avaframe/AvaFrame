#! /usr/bin/python

from aimectask import *

import re
import aimecdata
import aimecresult
import esrigrid
import numpy
import math

from PyQt5 import QtCore
from PyQt5 import QtGui

from coordtrans_process.processData import rasterize

class AimecTaskRunout(AimecTask):

    def __init__self(self):

        super(AimecTaskRunout, self).__init__()

    def cliName(self):

        return '-runout'

    def name(self):

        return 'Global Runout Comparison'

    def resultPath(self):

        return 'runout'

    def description(self):

        return ('Compare simulation outcomes with dokumented runout areas.<br/>'
               'Requires a <b>*runout.xyz</b> file with the dokuemnted runout area<br/>'
               'and a <b>*interest.xyz</b> with the area of intererest in the <b>doku input path</b> <br/>'
               'and pressure data files <b>*Id.txt</b> (from simulations) in the <b>pressure path</b>. <br/>'
               'Results will be written as *.txt')

    def validateData(self, data):

        #TODO fixme
        return True

    def run(self, data, callBack = None):

        #read measruemntresults
        dokuFileName = data.pathDocRunOut
        areaFileName = data.pathAOI

        with open(dokuFileName) as csvfile:
            dokudata = numpy.genfromtxt(csvfile, dtype = float)
        with open(areaFileName) as csvfile:
            areadata = numpy.genfromtxt(csvfile, dtype = float)
        dokuX = [dokudata[i,0] for i in range(len(dokudata))]
        dokuY = [dokudata[i,1] for i in range(len(dokudata))]
        areaX = [areadata[i,0] for i in range(len(areadata))]
        areaY = [areadata[i,1] for i in range(len(areadata))]
        areaPolygon = QtGui.QPolygonF([QtCore.QPointF(areaX[i], areaY[i]) for i in range(len(areaX))])
        dokuPolygon = QtGui.QPolygonF([QtCore.QPointF(dokuX[i], dokuY[i]) for i in range(len(dokuX))])
        boundingRect = areaPolygon.boundingRect()

        QtCore.qDebug('using doku %s/%s'%(dokuX, dokuY))
        QtCore.qDebug('using area of interest %s/%s'%(areaX, areaY))
        QtCore.qDebug('boundaryrect of interest %s'%(boundingRect))

        plim = data.calcPressureLimit
        QtCore.qWarning('assuming that pressure limt is %s '%(plim))
        rasterW = 5
        QtCore.qWarning('assuming that rasterwidth is %s '%(rasterW))

        # mapping for pixlcoords to real world coords
        globalTransform = QtGui.QTransform()
        globalTransform.translate(boundingRect.left(), boundingRect.top())
        globalTransform.scale(rasterW, rasterW)
        globalTransform = globalTransform.inverted()[0]

        pixlRect = globalTransform.mapRect(boundingRect)

        pixlRect.setHeight(math.ceil(pixlRect.height()))
        pixlRect.setWidth(math.ceil(pixlRect.width()))

        areaRaster = rasterize.rasterizePolygon(globalTransform.map(areaPolygon), pixlRect.size().toSize()).T
        dokuRaster = rasterize.rasterizePolygon(globalTransform.map(dokuPolygon), pixlRect.size().toSize()).T


        areaC = numpy.count_nonzero(areaRaster)
        doku = numpy.logical_and(areaRaster, dokuRaster)

        fillC = numpy.count_nonzero(doku)
        notfillC = areaC-fillC

        #create header
        #read simulationoutcomes
        result = aimecresult.AimecResult(self, aimecresult.AimecReference(
        [fillC, notfillC, fillC+notfillC], ['full', 'empty', 'sum'], [1, 2, 3]))

        fileList = QtCore.QDir(data.pathPressure).entryList(QtCore.QDir.Files)
        i = 0
        for fileName in fileList:
            if callBack:
                if callBack.getAbort():
                    QtCore.qWarning('RunOut task aborted by user.')
                    return None
            try:
                simulationId = int(re.search( r'(\d+)', fileName, re.M|re.I).group(1))
            except:
                simulationId = -1
            pressureGrid = esrigrid.esriOpen(data.pathPressure + QtCore.QString(QtCore.QDir.separator()) + fileName)

            simRaster = numpy.zeros(shape=(pixlRect.size().width(), pixlRect.size().height()))
            x_begin = int(math.floor((boundingRect.left()-pressureGrid.left())/rasterW))
            y_begin = int(math.floor((boundingRect.top()-pressureGrid.top())/rasterW))
            x_end = int(math.floor((boundingRect.right()-pressureGrid.left())/rasterW))
            y_end = int(math.floor((boundingRect.bottom()-pressureGrid.top())/rasterW))

            for x_index in range(x_end-x_begin):
                for y_index in range(y_end-y_begin):
                    if pressureGrid.valueAtXYCell(x_begin+x_index+1,  y_begin+y_index+1) > plim:
                        simRaster[x_index, y_index] = 1


            pos = numpy.logical_and(dokuRaster, simRaster)
            neg = numpy.logical_and(numpy.logical_not(dokuRaster), numpy.logical_not(simRaster))
            false_pos = numpy.logical_and(numpy.logical_not(dokuRaster), simRaster)
            false_neg = numpy.logical_and(dokuRaster, numpy.logical_not(simRaster))

            pos = numpy.count_nonzero(numpy.logical_and(pos, areaRaster))
            neg = numpy.count_nonzero(numpy.logical_and(neg, areaRaster))
            false_pos = numpy.count_nonzero(numpy.logical_and(false_pos, areaRaster))
            false_neg = numpy.count_nonzero(numpy.logical_and(false_neg, areaRaster))

            print(i, '->', pos, neg, false_pos, false_neg, pos + neg, false_pos + false_neg, pos + neg + false_pos + false_neg)

            result.simulations.append(aimecresult.AimecResultSingle(simulationId, [pos, neg, pos+neg], [false_pos, false_neg, false_pos+false_neg]))
            i += 1
            if callBack:
                callBack.setProgress(i, len(fileList))

        result.save(data.pathResult, self.resultPath())
        return result



if __name__ == '__main__':

    indata = aimecdata.AimecData()
    if False:
        indata.pathDocRunOut = '/home/matti/Frei/cosica/doku/vdls/vdls_103_doku_matti_runout.xyz'
        indata.pathPressure = '/home/matti/Frei/cosica/auto_results/samos/vdls_new_ss/pressure'
        indata.pathAOI = '/home/matti/Frei/cosica/doku/vdls/area_of_interest.xyz'
    else:
        indata.pathDocRunOut = '/home/matti/Frei/cosica/doku/rgf/rgf19970417_impacted_area_runout.xyz'
        indata.pathPressure = '/home/matti/Frei/cosica/auto_results/samos/rgf_new_ss/pressure'
        indata.pathAOI = '/home/matti/Frei/cosica/doku/rgf/area_of_interest.xyz'
    indata.pathResult = './testresults'
    indata.calcPressureLimit = 1.

    task = AimecTaskRunout()

    print('data is ok = ' + task.validateData(indata).__str__())

    print(task.run(indata).__str__())

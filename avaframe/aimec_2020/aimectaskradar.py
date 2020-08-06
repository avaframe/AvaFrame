#! /usr/bin/python

from aimectask import *

import re
import aimecdata
import aimecresult
import esrigrid
import numpy
import math

from PyQt5 import QtCore

class AimecTaskRadar(AimecTask):

    def __init__self(self):

        super(AimecTaskRadar, self).__init__()

    def cliName(self):

        return '-radar'

    def name(self):

        return 'Radar Comparison'

    def resultPath(self):

        return 'radar'

    def description(self):

        return ('Compare simulation outcomes with measurments from doppler radars.<br/>'
               'Requires a <b>*radarpath.xyz</b> file and a <b>*velocity.txt</b> file in the <b>radar input path</b> <br/>'
               'and pressure data files <b>*Id.txt</b> (from simulations) in the <b>pressure path</b>. <br/>'
               'Results will be written as *.txt')

    def validateData(self, data):

        #TODO fixme
        return True

    def run(self, data, callBack = None):

        #read measruemntresults
        velocityFileName = None
        pathFileName = None
        for fileName in QtCore.QDir(data.pathDocRadar).entryList(QtCore.QDir.Files):
            if 'velocity' in fileName:
                velocityFileName = data.pathDocRadar + QtCore.QString(QtCore.QDir.separator()) + fileName
            if 'path' in fileName:
                pathFileName = data.pathDocRadar + QtCore.QString(QtCore.QDir.separator()) + fileName
        if velocityFileName == None:
            raise RuntimeError('unable to find velocity file in %s'%(data.pathDocRadar))
        if pathFileName == None:
            raise RuntimeError('unable to find path file')

        with open(velocityFileName) as csvfile:
            velodata = numpy.genfromtxt(csvfile, dtype = float)
            radar = RadarMeasurment(velodata[:,1], velodata[:,3], velodata[:,5], velodata[:,4])
        with open(pathFileName) as csvfile:
             radarpos = numpy.genfromtxt(csvfile, dtype = float)
             radar.center = numpy.array([radarpos[0, 0], radarpos[0, 1]])
             radar.direction = numpy.array([radarpos[1,0], radarpos[1, 1]])-radar.center
        QtCore.qDebug('using range gates %s'%numpy.around(radar.x))
        QtCore.qDebug('using velocity %s'%numpy.around(radar.u))
        QtCore.qDebug('using radar position %s'%(numpy.around(radar.center)))
        QtCore.qDebug('using radar direction %s'%(numpy.around(radar.direction)))

        QtCore.qWarning('assuming that angle of radar is +-%s grad'%(radar.radarAngle))
        rho = data.density
        QtCore.qWarning('assuming that density is %s'%(rho))
        prefix = 1000
        QtCore.qWarning('assuming that pressure prefix is %s'%(prefix))

        #create header
        header = []
        i = 0
        for x in radar.x:
            header.append('u%s'%(i))
            i += 1

        result = aimecresult.AimecResult(self, aimecresult.AimecReference(
        radar.u, header, radar.x, radar.u_lower, radar.u_upper))
        result.reference.detailPos = numpy.zeros(shape=(int(numpy.linalg.norm(radar.direction)))/10)
        result.reference.detail = numpy.zeros(shape=(int(numpy.linalg.norm(radar.direction)))/10)
        for j in range(len(result.reference.detail)):
            result.reference.detail[j] = float('NAN')
            result.reference.detailPos[j] = 5+10*j
        #read simulationoutcomes
        cellRows = []
        cellCols = []
        cellPath = []
        rangeGates = []
        fileList = QtCore.QDir(data.pathPressure).entryList(QtCore.QDir.Files)
        i = 0
        for fileName in fileList:
            if callBack:
                if callBack.getAbort():
                    QtCore.qWarning('Radar task aborted by user.')
                    return None
            try:
                simulationId = int(re.search( r'(\d+)', fileName, re.M|re.I).group(1))
            except:
                simulationId = -1
            pressureGrid = esrigrid.esriOpen(data.pathPressure + QtCore.QString(QtCore.QDir.separator()) + fileName)
            maxVelocity = numpy.zeros(shape=(len(radar.x)))
            maxVelocityD = numpy.zeros(shape=(int(numpy.linalg.norm(radar.direction)))/10)
            #calculation of maximum velocity within a rangegate
            if len(cellRows) == 0:
                QtCore.qDebug('initializing mesh')
                for cell in pressureGrid:
                    if not radar.inSight(cell.center()):
                        continue
                    cellRows.append(cell.row)
                    cellCols.append(cell.col)
                    radarVec = cell.center()-radar.center
                    radarrange = numpy.linalg.norm(radarVec)
                    rangeGateId = radar.getRangeGate(radarrange)
                    cellPath.append(radarrange)
                    rangeGates.append(rangeGateId)
            for cellIndex in range(len(cellRows)):
                cell = esrigrid.EsriGridCell(pressureGrid, cellCols[cellIndex], cellRows[cellIndex])
                maxVelocity[rangeGates[cellIndex]] = max(maxVelocity[rangeGates[cellIndex]], cell.value())
                if not int(cellPath[cellIndex])/10 >= len(maxVelocityD):
                    maxVelocityD[int(cellPath[cellIndex])/10] = max(maxVelocityD[int(cellPath[cellIndex])/10], cell.value())
            maxVelocity = numpy.sqrt(maxVelocity*prefix/rho)
            maxVelocityD = numpy.sqrt(maxVelocityD*prefix/rho)
            error = abs(maxVelocity-radar.u)
            result.simulations.append(aimecresult.AimecResultSingle(simulationId, maxVelocity, error, maxVelocityD))
            i += 1
            if callBack:
                callBack.setProgress(i, len(fileList))
        result.save(data.pathResult, self.resultPath())
        return result

class RadarMeasurment(object):

    def __init__(self, x, u, u_lower, u_upper, pos = numpy.array([0, 0]), end = numpy.array([0,0])):
        self.x = x
        self.u = u
        self.u_lower = u_lower
        self.u_upper = u_upper
        self.center = pos
        self.direction = end-pos
        self.radarAngle = 15
        #calculating ranges of rangegates, lol
        self.ranges = []
        self.ranges.append([x[0]-(x[1]-x[0])/2., x[0]+(x[1]-x[0])/2.])
        for i in range(1,len(x)-1):
            self.ranges.append([x[i]-(x[i]-x[i-1])/2., x[i]+(x[i]-x[i-1])/2.])
        self.ranges.append([x[-1]-(x[-1]-x[-2])/2., x[-1]+(x[-1]-x[-2])/2.])

    def getRangeGate(self, x):
        i = 0
        if self.ranges[0][0] > x:
            return -1
        if self.ranges[-1][-1] < x:
            return -1
        for r in self.ranges:
            if r[0]< x and r[1] > x:
                return i
            i += 1
        return -1

    def inSight(self, cellCenter):
        radarVec = cellCenter-self.center
        a1 = numpy.dot(radarVec,self.direction)
        a2 = numpy.linalg.norm(numpy.cross(radarVec,self.direction))
        angle = numpy.arctan2(a2, a1)*180/math.pi
        if angle < self.radarAngle and angle > -self.radarAngle:
            return True
        return False

if __name__ == '__main__':


    a = numpy.array([1, 0])
    b = numpy.array([1, 1])
    angle1 = numpy.dot(a,b)
    angle2 = numpy.linalg.norm(numpy.cross(a,b))
    angle = numpy.arctan2(angle2, angle1)

    print(angle)

    indata = aimecdata.AimecData()
    indata.pathDocRadar = './radartestdata/measurment'
    indata.pathPressure = './radartestdata/pressure'
    indata.pathResult = './testresults'
    indata.density = 200.

    task = AimecTaskRadar()

    print('data is ok = ' + task.validateData(indata).__str__())

    print(task.run(indata).__str__())

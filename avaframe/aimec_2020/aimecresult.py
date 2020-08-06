#! /usr/bin/python

import math
import os
# from PyQt5 import QtCore

class AimecResultSingle(object):

    def __init__(self, simulationId, resultData, error, detail = None):

        self.simulationId = simulationId
        self.resultData = resultData
        self.error = error
        self.detail = detail

    def getError(self):

        errsum = 0
        for err in self.error:
            errsum += err**2
        errsum = errsum/len(self.error)
        return math.sqrt(errsum)

class AimecReference(object):

   def __init__(self, data, dataHeader, dataPos = None, dataMin = None, dataMax = None, detail = None, detailPos = None):

       self.data = data
       self.dataHeader = dataHeader
       self.dataPos = dataPos
       self.dataMin = dataMin
       self.dataMax = dataMax
       self.detail = detail
       self.detailPos = detailPos

class AimecResult(object):

    def __init__(self, task, reference):

        self.task = task
        self.simulations = []
        self.reference = reference

    def updateErrorsAuto(self):

        for sim in self.simulations:
            sim.error = [math.fabs(sim.resultData[i]-self.reference.data[i]) for i in range(len(self.reference.data))]


    def getNormErrors(self):
        """
        get an Error in form of a float:
        0 is the best simulation, 1 is the worst simulation
        """

        minErrorId, maxErrorId, minError, maxError = self.bestWorstSimulations()
        normfunction = lambda err: (err-minError)/(maxError-minError)
        return [normfunction(sim.getError()) for sim in self.simulations]

    def errorRanking(self):

        errorL = [sim.getError() for sim in self.simulations]
        return sorted(range(len(errorL)), key=errorL.__getitem__)

    def bestWorstSimulations(self):

        minError = float('inf')
        maxError = 0
        minErrorId = -1
        maxErrorId = -1
        i = 0
        for sim in self.simulations:
            error = sim.getError()
            if error < minError:
                minError = error
                minErrorId = i
            if error > maxError:
                maxError = error
                maxErrorId = i
            i += 1
        return minErrorId, maxErrorId, minError, maxError

    def bestSimulation(self):

        return self.bestWorstSimulations()[0]

    def __str__(self):

        string = 'results of %s - '%(self.task.name())
        bestIndex = self.bestSimulation()
        string += 'simulation %s (of %s) is the best (%s)'%(
            self.simulations[bestIndex].simulationId,
            len(self.simulations),
            self.simulations[bestIndex].getError())
        return string

    def save(self, pathname, taskpath):
        pathname = pathname+'/' + taskpath
        # path = QtCore.QDir(pathname)
        if not os.path.isdir(pathname):
            os.path.mkdir(pathname)
        else:
            print('[AimecResults] folder %s already exisits. Overwriting with new results.' %pathname)
        for mode in ['result', 'error', 'detail']:
            if 'detail' in mode:
                if not self.reference.detailPos is not None:
                    continue
            pathnameF = pathname+'/' + '%s.txt'%mode
            outfile = open(pathnameF, 'w')
            outfile.write('#standard aimec %s output from task "%s"\n' %(mode, self.task.name()))
            outfile.write('#line 1 (after comment) contains the header\n')
            outfile.write('#line 2 (after comment) contains the position of data\n')
            outfile.write('#line 3 (after comment) contains the referece data\n')
            outfile.write('#the following lines contain the simulationId and\n')
            outfile.write('#the calculated results from simulations\n')
            if 'error' in mode:
                outfile.write('#the last column contains the normed error\n')
            outfile.write('#number of simulations = %s\n'%(len(self.simulations)))
            outfile.write('#number of values to compare = %s\n'%(len(self.reference.data)))

            if 'detail' in mode:
                outfile.write('id;'+';'.join(['val %s'%i for i in range(len(self.reference.detail))])+'\n')
            elif 'error' in mode:
                outfile.write('id;'+';'.join(self.reference.dataHeader)+';err\n')
            else:
                outfile.write('id;'+';'.join(self.reference.dataHeader)+'\n')

            if 'detail' in mode:
                 #position of data
                outfile.write('0')
                for val in self.reference.detailPos:
                    outfile.write(';%8.4f' %val)
                outfile.write('\n')
                #reference
                outfile.write('0;')
                for val in self.reference.detail:
                    outfile.write(';%8.4f' %val)
                outfile.write('\n')

                #simulations
                for sim in self.simulations:
                    outfile.write('%s' %sim.simulationId)
                    for val in sim.detail:
                        outfile.write(';%8.4f' %val)
                    outfile.write('\n')
                outfile.close()

            else:
                #position of data
                outfile.write('0')
                for val in self.reference.dataPos:
                    outfile.write(';%8.4f' %val)
                if 'error' in mode:
                    outfile.write(';%8.4f' %float('NAN'))
                outfile.write('\n')
                #reference
                outfile.write('0')
                for val in self.reference.data:
                    outfile.write(';%8.4f' %val)
                if 'error' in mode:
                    outfile.write(';%8.4f' %float('NAN'))
                outfile.write('\n')

                #simulations
                for sim in self.simulations:
                    outfile.write('%s' %sim.simulationId)
                    if 'error' in mode:
                        for val in sim.error:
                            outfile.write(';%8.4f' %val)
                        outfile.write(';%8.4f\n'%sim.getError())
                    else:
                        for val in sim.resultData:
                            outfile.write(';%8.4f' %val)
                        outfile.write('\n')
                outfile.close()

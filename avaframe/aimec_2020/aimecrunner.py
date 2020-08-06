# -*- coding: utf-8 -*-

import aimectask
import aimectaskradar
import aimectaskrunout
import aimectaskcoordtrans
import copy

class AimecRunner(object):


    # def __init__(self, parent = None):
    def __init__(self):
        # super(AimecRunner, self).__init__(parent)

        self.tasks = [aimectaskcoordtrans.AimecTaskCoordTrans(),
                      aimectaskradar.AimecTaskRadar(),
                      aimectaskrunout.AimecTaskRunout()]

    def __del__(self):
        pass

    def setProgress(self, progress, overall):
        print('[ARunner] progressChanged',progress)

    def run(self):

        [data, tasklist] = self.toDo
        self.toDo = 0

        self.results = []
        for task in self.tasks:
            print('[ARunner] taskChanged',task)

            if task in tasklist:
                print('[ARunner] validating %s'%(task.name()))
                if task.validateData(data):
                    print('[ARunner] running %s'%(task.name()))
                    try:
                        self.results.append(task.run(data, self))
                        print('[ARunner] finished %s: %s'%(task.name(), self.results[-1]))
                    except RuntimeError as e:
                        print('[ARunner] %s raised a runtime exception: %s'%(task.name(), e.message))
                        self.results.append(None)
                else:
                    print('[ARunner] validation of %s failed. Skipping task'%(task.name()))
                    self.results.append(None)
            else:
                self.results.append(None)
        print('[ARunner] finished')
        self.results = []
        return

    def getAbort(self):
        return False

    def abortCalculation(self):
        pass

    def runall(self, data, taskList):
        self.toDo = [copy.copy(data), copy.copy(taskList)]
        self.run()

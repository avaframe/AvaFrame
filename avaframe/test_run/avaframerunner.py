import avaframetaskinit
import avaframetaskalphabeta
import copy
import logging
logmain = logging.getLogger(__name__)

class AvaframeRunner(object):

    # def __init__(self, parent = None):
    def __init__(self):
        # super(AvaframeRunner, self).__init__(parent)

        self.tasks = [avaframetaskinit.AvaframeTaskInitialize(),
                      avaframetaskalphabeta.AvaframeTaskAlphaBeta()]

    def __del__(self):
        pass

    def setProgress(self, progress, overall):
        logmain.info('progressChanged', progress)

    def run(self):

        [data, tasklist] = self.toDo
        self.toDo = 0

        self.results = []
        for task in self.tasks:


            if task in tasklist:
                logmain.info('taskChanged %s' % (task.name()))
                logmain.info('validating %s' % (task.name()))
                if task.validateData(data):
                    logmain.info('running %s' % (task.name()))
                    try:
                        self.results.append(task.run(data, self))
                        # fh = logging.FileHandler("Main_Avaframe.log")
                        # logmain.addHandler(fh)
                        logmain.info('finished %s' % (task.name()))
                    except RuntimeError as e:
                        logmain.info('%s raised a runtime exception: %s' %
                              (task.name(), e.message))
                        self.results.append(None)
                else:
                    logmain.info('validation of %s failed. Skipping task' % (task.name()))
                    self.results.append(None)
            # else:
            #     logmain.info('%s is not a task' % (task.name()))
            #     self.results.append(None)
        logmain.info('finished')
        self.results = []
        return

    def getAbort(self):
        return False

    def abortCalculation(self):
        pass

    def runall(self, data, taskList):
        self.toDo = [copy.copy(data), copy.copy(taskList)]
        self.run()

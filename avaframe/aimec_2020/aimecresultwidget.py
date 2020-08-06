from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

from PyQt4 import QtCore
from PyQt4 import QtGui


import aimecresult

class AimecResultWidget(QtGui.QWidget):

    def __init__(self, result, parent=None):

        super(AimecResultWidget, self).__init__(parent)

        self.result = result

        self.figure = Figure()
        self.canvas = FigureCanvasQTAgg(self.figure)

        self.axis = self.figure.add_subplot(111)

        self.layout = QtGui.QGridLayout(self)

        self.plotSelector = QtGui.QComboBox(self)
        self.plotSelector.addItem('points')
        self.plotSelector.addItem('simulations')

        self.itemSelector = QtGui.QComboBox(self)
        self.referenceSelector = QtGui.QCheckBox('reference', self)
        self.referenceSelector.setChecked(True)
        self.detailSelector = QtGui.QCheckBox('detailed', self)
        self.detailSelector.setChecked(True)

        self.layout.addWidget(QtGui.QLabel('show on x-axis:', self), 0, 0, 1, 1)
        self.layout.addWidget(self.plotSelector, 0, 1, 1, 1)
        self.layout.addWidget(QtGui.QLabel('show simulation/position:', self), 0, 2, 1, 1)
        self.layout.addWidget(self.itemSelector, 0, 3, 1, 1)
        self.layout.addWidget(self.referenceSelector, 0, 4, 1, 1)
        self.layout.addWidget(self.detailSelector, 2, 4, 1, 1)
        self.layout.addWidget(self.canvas, 1, 0, 20, 4)

        for label in self.findChildren(QtGui.QLabel):
            label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)

        self.plotInit(0)

        QtCore.QObject.connect(self.plotSelector, QtCore.SIGNAL('currentIndexChanged(int)'), self.plotInit)
        QtCore.QObject.connect(self.itemSelector, QtCore.SIGNAL('currentIndexChanged(int)'), self.updatePlot)

        QtCore.QObject.connect(self.referenceSelector, QtCore.SIGNAL('stateChanged(int)'), self.updatePlot)
        QtCore.QObject.connect(self.detailSelector, QtCore.SIGNAL('stateChanged(int)'), self.updatePlot)

    def updatePlot(self):

        self.axis.clear()

        mode = self.plotSelector.currentIndex()
        showReference = self.referenceSelector.isChecked()
        item = self.itemSelector.currentIndex()
        detailed = self.detailSelector.isChecked()
        if mode == 0: #points on x
            if showReference:
                if self.result.reference.dataMax is not None:
#                     self.axis.errorbar(self.result.reference.dataPos, self.result.reference.data,
#                                        yerr = [self.result.reference.dataMin,self.result.reference.dataMax])
                     self.axis.errorbar(self.result.reference.dataPos, self.result.reference.data,
                                        yerr = [[self.result.reference.data[i]-self.result.reference.dataMin[i] for i in range(len(self.result.reference.data))],
                                                [self.result.reference.dataMax[i]-self.result.reference.data[i] for i in range(len(self.result.reference.data))]])
                else:
                    self.axis.plot(self.result.reference.dataPos, self.result.reference.data)
            if self.itemSelector.currentIndex() == 0: #onlyreference
                plotRange = []
            elif self.itemSelector.currentIndex() == 1:#all
                plotRange = range(len(self.result.simulations))
            elif self.itemSelector.currentIndex() == 2:#best
                plotRange = [self.result.errorRanking()[0]]
            else:
                plotRange = self.result.errorRanking()[0:5]
            for index in plotRange:
                if detailed:
                    self.axis.plot(self.result.reference.detailPos, self.result.simulations[index].detail)
                else:
                    self.axis.plot(self.result.reference.dataPos, self.result.simulations[index].resultData)
        elif mode == 1: # simulations on x
            if showReference:
                self.axis.plot([self.result.simulations[i].simulationId for i in range(len(self.result.simulations))],
                                [self.result.reference.data[item] for i in range(len(self.result.simulations))])
            self.axis.plot([self.result.simulations[i].simulationId for i in range(len(self.result.simulations))],
                           [self.result.simulations[i].resultData[item] for i in range(len(self.result.simulations))])

        self.canvas.draw()

    def plotInit(self, mode):

        if mode == 0:
            self.itemSelector.clear()
            self.itemSelector.addItem('only reference')
            self.itemSelector.addItem('all')
            self.itemSelector.addItem('best')
            self.itemSelector.addItem('best5')
        else:
            self.itemSelector.clear()
            for header in self.result.reference.dataHeader:
                self.itemSelector.addItem('%s'%header)

        self.updatePlot()

if __name__ == "__main__":
    import sys
    import aimectaskradar

    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('MyWindow')

    r = aimecresult.AimecResult(aimectaskradar.AimecTaskRadar(),
                                aimecresult.AimecReference([5, 10, 20, 10],
                                                           ['u1', 'u2', 'u3', 'u4'],
                                                           [1, 2, 3, 4], [2, 7, 15, 3], [12, 13, 23, 15]))

    r.simulations.append(aimecresult.AimecResultSingle(1, [2, 3, 4, 5], []))
    r.simulations.append(aimecresult.AimecResultSingle(2, [2, 3, 6, 5], []))
    r.simulations.append(aimecresult.AimecResultSingle(3, [3, 4, 4, 5], []))
    r.updateErrorsAuto()

    main = AimecResultWidget(r)
    main.show()

    sys.exit(app.exec_())

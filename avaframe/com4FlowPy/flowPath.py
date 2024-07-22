import numpy as np
import math
import csv

class Path:
    '''Class contains a path, containing one startcell  and corresponding child cells'''

    def __init__(self, dem, startcellRow, startcellCol, genList):
        self.dem = dem
        self.cellsize = genList[0][0].cellsize
        self.alpha = genList[0][0].alpha
        self.exp = genList[0][0].exp
        self.maxZDelta = genList[0][0].max_z_delta
        self.genList = genList
        self.startcellRow = startcellRow
        self.startcellCol = startcellCol

        self.dropHeight = 0
        self.travelLength = 0

        self.zDeltaArray = np.zeros_like(self.dem, dtype=np.float32)
        self.flowEnergyArray = np.zeros_like(self.dem, dtype=np.float32)
        self.fluxArray = np.zeros_like(self.dem, dtype=np.float32)
        '''
        self.travel_length_array = np.zeros_like(self.dem, dtype=np.float32)
        self.generation_array = np.full_like(self.dem, np.nan, dtype=np.float32)
        '''

        self.zDeltaGeneration = []
        self.fluxGeneration = []
        self.travelLengthGeneration = []
        self.flowEnergyGeneration = []
        self.rowGeneration = []
        self.colGeneration = []
        self.altitudeGeneration = []
        self.gammaGeneration = []
        #self.pathArea = 0
        self.flux_gen = []

    
    def getVariablesGeneration(self):
        '''
            write LISTS with size and format of genList 
            (the main list contains lists for every generation)
        '''
        for gen, cellList in enumerate(self.genList):
            cellListZDelta = []
            cellListFlux = []
            cellListMinDistance = []
            cellListFlowEnergy = []
            cellListRow = []
            cellListCol = []
            cellListAlt = []
            cellListGamma = []
            cellListFlux_gen = []

            for cell in cellList:
                cellListZDelta.append(cell.z_delta)
                cellListFlux.append(cell.flux)
                cellListMinDistance.append(cell.min_distance)
                cellListFlowEnergy.append(cell.flowEnergy)
                cellListRow.append(cell.rowindex)
                cellListCol.append(cell.colindex)
                cellListAlt.append(cell.altitude)
                cellListGamma.append(cell.max_gamma)
                #cellListFlux_gen.append(cell.flux_generation)

            self.zDeltaGeneration.append(cellListZDelta)
            self.fluxGeneration.append(cellListFlux)
            self.travelLengthGeneration.append(cellListMinDistance)
            self.flowEnergyGeneration.append(cellListFlowEnergy)
            self.rowGeneration.append(cellListRow)
            self.colGeneration.append(cellListCol)
            self.altitudeGeneration.append(cellListAlt)
            self.gammaGeneration.append(cellListGamma)
            #self.flux_gen.append(cellListFlux_gen)
        
    
    def getPathArrays(self):
        ''' 
            write ARRAYS with size of dem containing the variable values of every path
            value 0 means, the path does not hit the cell
        '''
        for gen, cellList in enumerate(self.genList):
            for cell in cellList:
                self.zDeltaArray[cell.rowindex, cell.colindex] = max(self.zDeltaArray[cell.rowindex, cell.colindex], cell.z_delta)
                self.flowEnergyArray[cell.rowindex, cell.colindex] = max(self.flowEnergyArray[cell.rowindex, cell.colindex], cell.flowEnergy)
                self.fluxArray[cell.rowindex, cell.colindex] = max(self.fluxArray[cell.rowindex, cell.colindex], cell.flux)
                
                '''
                self.travel_length_array[cell.rowindex, cell.colindex] = max(self.travel_length_array[cell.rowindex, cell.colindex], cell.min_distance)
                self.generation_array[cell.rowindex, cell.colindex] = gen
                '''

    
    def calcThalwegCenterof(self, variable, variableCo):
        '''
        Parameters:

        variable: list
            variable, which is centered (in format genList)
        variableCo: list
            center of variableCo is calculated (variable is weighted) (in format genList)

        Return:

        variableSum: np.array
            sum of variable per generation
        coVar: np.array
            centered variable (per generation)
        '''
        coVar = np.zeros(len(self.genList))
        variableSum = np.zeros(len(self.genList))
        for gen in range(0,len(self.genList)):
            var = np.array(variable[gen])
            co = np.array(variableCo[gen])
            variableSum[gen] = np.sum(var)
            variableCoSum = np.sum(co)
            if variableCoSum > 0: # flow_energy is 0 in generation 0
                coVar[gen] = 1 / variableCoSum * np.sum(var * co)
            else:
                coVar[gen] = np.sum(var)
        return variableSum, coVar


    def getCenterofs(self, variables):
        '''
        calculate sum of variable for every iteration step/ generation and
        center of energy & flux for the following variables:

        Parameters:

        variables: list
            List of variables that should be weighted (with center of energy and flux)
        '''
        self.getVariablesGeneration()

        for varName in variables:
            values = getattr(self, f'{varName}Generation')
            sumF, coF = self.calcThalwegCenterof(values, self.fluxGeneration) # center of flux of every variable
            sumE, coE = self.calcThalwegCenterof(values, self.flowEnergyGeneration) # center of energy of every variable
            sumZd, coZd = self.calcThalwegCenterof(values, self.zDeltaGeneration) # center of energy of every variable


            setattr(self, f'{varName}SumGen', sumF)
            setattr(self, f'{varName}CoF', coF)
            setattr(self, f'{varName}CoE', coE)
            setattr(self, f'{varName}CoZd', coZd)


    def calcRunoutAngle(self):
        '''
        Calculate runout angle of thalweg
        '''

        self.dropHeightCoE = max(self.altitudeCoE)-min(self.altitudeCoE)
        self.travelLengthCoE = max(self.sCoE)
        if self.dropHeightCoE > 0:
            self.runoutAngleCoE = np.rad2deg(np.arctan(self.dropHeightCoE / self.travelLengthCoE))

        self.dropHeightCoF = max(self.altitudeCoF)-min(self.altitudeCoF)
        self.travelLengthCoF = max(self.sCoF)
        if self.dropHeightCoF > 0:
            self.runoutAngleCoF = np.rad2deg(np.arctan(self.dropHeightCoF / self.travelLengthCoF))


    def saveThalwegData(self, variable, variableValue, saveDir):
        '''
        Parameters:

        variable: str
            name of variable that is saved
        
        variableValue: np.array
            Array of thalweg data that is saved

        saveDir: pathlib.PosixPath
            directory, in which the thalweg data is saved

        '''

        with open(saveDir / (f"{variable}.csv"), mode='a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(variableValue)


    def saveAllThalwegData(self, saveDir):   
        variables = ['col', 'row', 'flux', 'flowEnergy', 'altitude', 'travelLength', 'zDelta', 'gamma']

        self.getCenterofs(variables)

        for varName in variables:
            value = getattr(self, f'{varName}CoE')
            self.saveThalwegData(varName, value, saveDir)

        metaData = {'alpha': round(self.alpha,1),
                    'exponent': self.exp,
                    'zDeltaMax': round(self.maxZDelta,1),
                    }
        np.save(saveDir / (f"meta.npy"), metaData)

        self.saveDict(saveDir)

        '''
        f = open(saveDir / (f"meta.csv"), mode='a')
        f.write(f'alpha: {round(self.alpha,2)}° \n \
                exponent: {self.exp}  \n \
                max zDelta: {round(self.maxZDelta,2)} m \n \
                v max: {round(np.sqrt(self.maxZDelta * 2 * 9.81),2)} m/s \n')
        f.close()
        '''

    def saveDict(self, saveDir):
        """
        save thalweg data. (One file per thalweg)
        Parameters:

        saveDir: pathlib.PosixPath
            directory, in which the thalweg data is saved

        """
        thalwegData = {'alpha': round(self.alpha,1),
                    'exponent': self.exp,
                    'zDeltaMax': round(self.maxZDelta,1),
                    }
        variables = ['col', 'row', 'flux', 'flowEnergy', 'altitude', 'travelLength', 'zDelta', 'gamma']
        centerOfs = ['CoE', 'CoF', 'CoZd']
        for varName in variables:
            for co in centerOfs:
                value = getattr(self, f'{varName}{co}')
                thalwegData[f'{varName}{co}'] = value

        self.getPathArrays()
        thalwegData['flowEnergyArray'] = self.flowEnergyArray
        thalwegData['zDeltaArray'] = self.zDeltaArray
        thalwegData['fluxArray'] = self.fluxArray

        self.stoppingCriteris()
        thalwegData['StoppingAlpha'] = self.StoppingAlpha
        thalwegData['StoppingVmax'] = self.StoppingVmax
        np.save(saveDir / (f"thalwegData_{self.startcellRow}_{self.startcellCol}.npy"), thalwegData)


def stoppingCriteria(self):
    self.calcRunoutAngle()
    if self.runoutAngleCoE == self.alpha:
        self.StoppingAlpha = True
    else:
        self.StoppingAlpha = False

    if max(self.zDeltaCoE) == self.maxZDelta:  #TODO: max(self.zDeltaCoE) oder max(self.zDeltaGeneration)??
        self.StopppingVmax = True
    else:
        self.StoppingVmax = False
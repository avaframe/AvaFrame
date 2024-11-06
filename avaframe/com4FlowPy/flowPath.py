import numpy as np
import math
import csv
import pickle

class Path:
    '''Class contains a path, containing one startcell  and corresponding child cells'''

    def __init__(self, dem, startcellRow, startcellCol, genList, rasterAttributes):
        """ initializes a GMF path, that belongs to a startcell

        Parameters
        ----------
        dem: numpy array
            Digital elevation model
        startcellRow: int
            Row index of startcell
        startcellCol: int
            Column index of startcell
        genList: list
            contains all cells that belong to the path (per generation an extra list)
        rasterAttributes: dict
            contains information about the input rasters
        """
        self.dem = dem
        self.cellsize = rasterAttributes["cellsize"]
        self.xllcorner = rasterAttributes["xllcenter"] - self.cellsize / 2
        self.yllcorner = rasterAttributes["yllcenter"] - self.cellsize / 2
        self.nrows = rasterAttributes["nrows"]
        self.crs = rasterAttributes["crs"]

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

    def indizesToCoords(self, cols, rows):
        """calculates the row and column indices to the x and y coordinates

        Parameters
        ----------
        cols: numpy array
            column indices of cells belonging to path
        rows: numpy array
            row indices of cells belonging to path

        Returns
        ----------
        x: numpy array
            x coordinates (in m) of cells belonging to path
        y: numpy array
            y coordinates (in m) of cells belonging to path
        """
        x = cols * self.cellsize + self.xllcorner
        y = self.yllcorner + (self.nrows - rows) * self.cellsize
        return x,y


    def getVariablesGeneration(self):
        '''write lists with size and format of genList containing specific parameters
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
        '''write arrays with size of dem containing the maximum of the variable values of every path
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
        '''calculates for a specific variable the center of a specific variable (thalweg)

        Parameters
        ----------
        variable: list
            variable, which is centered (in format genList)
        variableCo: list
            center of variableCo is calculated (variable is weighted) (in format genList)

        Returns
        ----------
        variableSum: numpy array
            sum of variable per generation
        coVar: numpy array
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
        center of energy, flux and zDelta for the following variables:

        Parameters
        ----------
        variables: list
            List of variables that should be weighted (with center of energy and flux)
        '''
              
        self.getVariablesGeneration()

        for varName in variables:
            if varName in ['s', 'z', 'x', 'y', 'flowEnergyArray', 'zDeltaArray', 'fluxArray']:
                continue
            values = getattr(self, f'{varName}Generation')
            sumF, coF = self.calcThalwegCenterof(values, self.fluxGeneration) # center of flux of every variable
            sumE, coE = self.calcThalwegCenterof(values, self.flowEnergyGeneration) # center of energy of every variable
            sumZd, coZd = self.calcThalwegCenterof(values, self.zDeltaGeneration) # center of energy of every variable

            setattr(self, f'{varName}SumGen', sumF)
            setattr(self, f'{varName}CoF', coF)
            setattr(self, f'{varName}CoE', coE)
            setattr(self, f'{varName}CoZd', coZd)


    def saveDict(self, saveDir, centerOfs, variables):
        """
        save thalweg data. (One file per thalweg)

        Parameters
        ------------
        saveDir: pathlib.PosixPath
            directory, in which the thalweg data is saved
        centerOfs: list
            contains the center-of-variable names that are saved
        variables: list
            contains the variable names that are saved
        """
        thalwegData = {'alpha': round(self.alpha,1),
                    'exponent': self.exp,
                    'zDeltaMax': round(self.maxZDelta,1),
                    'crs': self.crs,
                    }
        variables = variables
        centerOfs = centerOfs
        
        for co in centerOfs:
            for varName in variables:
                if varName in ['x', 'y']:
                    # compute x and y coordinates of thalweg
                    x, y = self.indizesToCoords(getattr(self, f'col{co}'), getattr(self, f'row{co}'))
                    setattr(self, f'x{co}', x)
                    setattr(self, f'y{co}', y) 

                if varName in ['flowEnergyArray', 'zDeltaArray', 'fluxArray']:
                    self.getPathArrays()
                    value = getattr(self, f'{varName}')
                elif varName == 'z':
                    value = getattr(self, f'altitude{co}')
                elif varName == 's':
                    value = getattr(self, f'travelLength{co}')
                else:
                    value = getattr(self, f'{varName}{co}')
                thalwegData[f'{varName}'] = value

            with open(saveDir / (f"thalwegData_{co}_{self.startcellRow}_{self.startcellCol}.pickle"), 'wb') as handle:
                pickle.dump(thalwegData, handle, protocol=pickle.HIGHEST_PROTOCOL)


    def calcAndSaveThalwegData(self, thalwegParameters):
        """main function for paths & thalwegs: calculates the thalweg and saves the data
        
        Parameters:
        ------------
        thalwegParameters: dict
            contains information to calculate and save the thalweg data (from .ini file)
        """
        saveDir = thalwegParameters["thalwegDir"]
        cos = eval(thalwegParameters["thalwegCenterOf"])
        variables = eval(thalwegParameters["thalwegVariables"])
        centerOfs = []
        for co in cos:
            co.lower()
            if co in ['energy', 'coe']:
                centerOf = 'CoE'
            elif co in ['flux', 'cof']:
                centerOf = 'CoF'
            elif co in ['zdelta', 'cozd']:
                centerOf = 'CoZd'
            centerOfs.append(centerOf)


        if 's' in variables:
            variables.append('travelLength')
        if 'z' in variables:
            variables.append('altitude')
        if 'x' in variables or 'y' in variables:
            variables.append('col')
            variables.append('row')            
        
        self.getCenterofs(variables)
        self.saveDict(saveDir, centerOfs, variables)
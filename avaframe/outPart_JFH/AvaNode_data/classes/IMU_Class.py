import numpy as np
import pandas as pd

class ImuData(object):

    def __init__(self):

        self.path = ''
        self.work_dir = ''
        self.header = ''
        self.gravi = 9.80665  #ToDO: make this variable depending on where we are (GPS)
        self.date = 0
        self.acc_factor = 0
        self.gyro_factor = 0
        self.mag_factor = 0
        self.data = 0
        
        self.is_synced = False
        self.sync_counter = []

    def read_data(self):
        print("Reading data!")
        leader_file = self.path
        header = ''
        mpuSettings = ""

        with open(leader_file, 'r') as file:
            for idx, line in enumerate(file):
                if idx < 5:
                    header = header + line[1:]

                if idx == 1:
                    print(line[1:-1])

                if idx == 3:
                    line = line[1:-2]  # get rid of comment markers
                    self.comment = line.split(',')

                if idx == 4:
                    line = line[1:-2]  # get rid of comment markers
                    mpuSettings = line.split(',')
           
            self.header = header
            self.date = mpuSettings[0]
            self.acc_factor = float(mpuSettings[1])
            self.gyro_factor = float(mpuSettings[2])
            self.mag_factor = float(mpuSettings[3])
            self.data = np.loadtxt(leader_file, comments='#', delimiter=',')            
            # define acc, gyro and mag vectors
            self.time = self.data[:, 0] / 1000000
            self.time = self.time - self.time[0]
            self.acc_X = self.data[:, 1] / self.acc_factor * self.gravi
            self.acc_Y = self.data[:, 2] / self.acc_factor * self.gravi
            self.acc_Z = self.data[:, 3] / self.acc_factor * self.gravi
            self.acc_tot = np.sqrt(self.acc_X ** 2 + self.acc_Y ** 2 + self.acc_Z ** 2)
            self.acc = np.column_stack((self.acc_X, self.acc_Y, self.acc_Z))
            self.gyro_X = self.data[:, 4] / self.gyro_factor
            self.gyro_Y = self.data[:, 5] / self.gyro_factor
            self.gyro_Z = self.data[:, 6] / self.gyro_factor
            self.gyro = np.column_stack((self.gyro_X, self.gyro_Y, self.gyro_Z))
            self.mag_X = self.data[:, 7] / self.mag_factor
            self.mag_Y = self.data[:, 8] / self.mag_factor
            self.mag_Z = self.data[:, 9] / self.mag_factor
            self.mag_tot = np.sqrt(self.mag_X ** 2 + self.mag_Y ** 2 + self.mag_Z ** 2)
            self.mag = np.column_stack((self.mag_X, self.mag_Y, self.mag_Z))
            self.sync = self.data[:, 10]
            # Datastructure for Calibrated values
            self.acc_cali = np.zeros_like(self.acc)
            self.acc_X_cali = np.zeros_like(self.acc_X)
            self.acc_Y_cali = np.zeros_like(self.acc_X)
            self.acc_Z_cali = np.zeros_like(self.acc_X)
            # Datastrucutre for Earth Frame values
            self.aEF = np.zeros_like(self.acc)
            self.axEF = np.zeros_like(self.acc_X)
            self.ayEF = np.zeros_like(self.acc_X)
            self.azEF = np.zeros_like(self.acc_X)
            self.atotEF = np.zeros_like(self.acc_X)
            self.vEF = np.zeros_like(self.acc)
            self.vxEF = np.zeros_like(self.acc_X)
            self.vyEF = np.zeros_like(self.acc_X)
            self.vzEF = np.zeros_like(self.acc_X)
            self.vtotEF = np.zeros_like(self.acc_X)
            self.pEF = np.zeros_like(self.acc)
            self.pxEF = np.zeros_like(self.acc_X)
            self.pyEF = np.zeros_like(self.acc_X)
            self.pzEF = np.zeros_like(self.acc_X)
            self.ptotEF = np.zeros_like(self.acc_X)

        print("Finished")
        
    def read_data_pd(self):
        """Reads File into a Pandas DataFrame."""
        print("Reading data!")
        leader_file = self.path
        header = ''
        mpuSettings = ""
        # Check for file version
        # Version 1 = Old, till 2022
        # Version 2 = New, from 2022 on
        
        with open(leader_file, 'r') as file:
            first_line = file.readline()
            header = first_line
            if not 'filename' in first_line and 'Data' in first_line: # Version 1
                version = 1
                print("Leader file version 1.")
            if 'filename' in first_line:
                version = 2
                print("Leader file version 2")
            if 'MicrocontrollerV2022.2' in first_line:
                version = 3
                print("Leader file version 3")
        
            if version == 1:
                for idx, line in enumerate(file):
                
                    if idx < 4:
                        header = header + line[1:]
    
                    if idx == 0:
                        print(line[1:-1])
    
                    if idx == 2:
                        line = line[1:-2]  # get rid of comment markers
                        self.comment = line.split(',')
    
                    if idx == 3:
                        line = line[1:-2]  # get rid of comment markers
                        mpuSettings = line.split(',')
                        
                    if idx > 4:
                        break

                self.header = header[:-1]
                self.date = mpuSettings[0]
                self.acc_factor = float(mpuSettings[1])
                self.gyro_factor = float(mpuSettings[2])
                self.mag_factor = float(mpuSettings[3])
                columns = ['time', 'acc_x', 'acc_y', 'acc_z', 'gyro_x', 'gyro_y', 'gyro_z', 'mag_x', 'mag_y', 'mag_z', 'sync']
                skiprows = 5
                
            elif version ==2:
                for idx, line in enumerate(file):
                
                    if idx < 12:
                        header = header + line[1:]
    
                    if idx == 2:
                        print(line[1:-1])
    
                    if idx == 3:
                        line = line[1:-2]  # get rid of comment markers
                        self.comment = line.split(',')
    
                    if idx == 11:
                        line = line[1:-2]  # get rid of comment markers
                        mpuSettings = line.split(',')

                self.header = header[:-1]
                self.date = mpuSettings[0]
                #self.acc_factor = float(mpuSettings[1])
                #self.gyro_factor = float(mpuSettings[12])
                #self.mag_factor = float(mpuSettings[15])
                self.acc_factor = 2048
                self.gyro_factor = 16.4
                self.mag_factor = 1/0.15
                columns = ['time', 'acc_x', 'acc_y', 'acc_z', 'gyro_x', 'gyro_y', 'gyro_z', 'mag_x', 'mag_y', 'mag_z', 'sync']
                skiprows = 13
                
            elif version ==3:
                for idx, line in enumerate(file):
                
                    if idx < 13:
                        header = header + line[1:]
    
                    if idx == 2:
                        print(line[1:-1])
    
                    if idx == 3:
                        line = line[1:-2]  # get rid of comment markers
                        self.comment = line.split(',')
    
                    if idx == 11:
                        line = line[1:-2]  # get rid of comment markers
                        mpuSettings = line.split(',')

                self.header = header[:-1]
                self.date = mpuSettings[0]
                #self.acc_factor = float(mpuSettings[1])
                #self.gyro_factor = float(mpuSettings[12])
                #self.mag_factor = float(mpuSettings[15])
                self.acc_factor = 2048
                self.gyro_factor = 16.4
                self.mag_factor = 1/0.15
                columns = ['time', 'acc_x', 'acc_y', 'acc_z', 'gyro_x', 'gyro_y', 'gyro_z', 'mag_x', 'mag_y', 'mag_z', 'sync']
                skiprows = 14
            
            self.data = pd.read_csv(leader_file, skiprows=skiprows, delimiter=',', names=columns, dtype={'time': int}, skipfooter=1, engine='python')
            self.data = self.data[:-1]
            self.data['time'] = self.data['time'].astype(int)
            # define acc, gyro and mag vectors
            for idx, time in enumerate(self.data['time']):
                if type(time) == str:
                    print("String on index {}".format(idx))
            self.time = self.data['time'].sub(self.data['time'][0])/1000000 # seconds
            #self.time = self.time - self.time[0]
            self.acc_X = self.data['acc_x'].astype(float) / self.acc_factor * self.gravi
            self.acc_Y = self.data['acc_y'].astype(float) / self.acc_factor * self.gravi
            self.acc_Z = self.data['acc_z'].astype(float) / self.acc_factor * self.gravi
            self.acc_tot = np.linalg.norm((self.acc_X, self.acc_Y, self.acc_Z), axis=0)
            self.acc = np.column_stack((self.acc_X, self.acc_Y, self.acc_Z))
            self.gyro_X = self.data['gyro_x'].astype(float) / self.gyro_factor
            self.gyro_Y = self.data['gyro_y'].astype(float) / self.gyro_factor
            self.gyro_Z = self.data['gyro_z'].astype(float) / self.gyro_factor
            self.gyro_tot = np.linalg.norm((self.gyro_X, self.gyro_Y, self.gyro_Z), axis=0)
            self.gyro = np.column_stack((self.gyro_X, self.gyro_Y, self.gyro_Z))
            self.mag_X = self.data['mag_x'].astype(float) / self.mag_factor
            self.mag_Y = self.data['mag_y'].astype(float) / self.mag_factor
            self.mag_Z = self.data['mag_z'].astype(float) / self.mag_factor
            self.mag_tot = np.sqrt(self.mag_X ** 2 + self.mag_Y ** 2 + self.mag_Z ** 2)
            self.mag = np.column_stack((self.mag_X.to_numpy(), self.mag_Y.to_numpy(), self.mag_Z.to_numpy()))
            self.sync = self.data['sync']
            
            for i in range(len(self.mag)):
                self.mag[i,:] = self.RotXYZ2RotationMatrix([180*np.pi/180,0,-90*np.pi/180])@self.mag[i,:]
                self.mag[i,:] = self.RotXYZ2RotationMatrix([0, 0, np.pi])@self.mag[i,:]
            
            self.mag_X = self.mag[:, 0]
            self.mag_Y = self.mag[:, 1]
            self.mag_Z = self.mag[:, 2]
            
        print("Finished")
    
    def RotationMatrix2RotXYZ(self, rotationMatrix):
        R=rotationMatrix
        #rot=np.array([0,0,0])
        rot=[0,0,0]
        rot[0] = np.arctan2(-R[1][2], R[2][2])
        rot[1] = np.arctan2(R[0][2], np.sqrt(abs(1. - R[0][2] * R[0][2]))) #fabs for safety, if small round up error in rotation matrix ...
        rot[2] = np.arctan2(-R[0][1], R[0][0])
        return np.array(rot)
    
    def RotXYZ2RotationMatrix(self, rot):
        c0 = np.cos(rot[0])
        s0 = np.sin(rot[0])
        c1 = np.cos(rot[1])
        s1 = np.sin(rot[1])
        c2 = np.cos(rot[2])
        s2 = np.sin(rot[2])
        
        return np.array([[c1*c2,-c1 * s2,s1],
                      [s0*s1*c2 + c0 * s2, -s0 * s1*s2 + c0 * c2,-s0 * c1],
                      [-c0 * s1*c2 + s0 * s2,c0*s1*s2 + s0 * c2,c0*c1 ]])
    
    def sync_leader(self):
        print("Syncing the Signals...")

        startvalue = self.sync[0]
        leader_peaks = 0
        leader_sync_counter = []

        counter_even = 0
        counter_odd = 1

        for value in self.sync:
            if value == startvalue:
                if value == 0:
                    leader_sync_counter.append(counter_even)
                if value == 1:
                    leader_sync_counter.append(counter_odd)
            if value != startvalue:
                leader_peaks += 1
                if value == 0:
                    counter_even += 2
                    leader_sync_counter.append(counter_even)
                if value == 1:
                    if leader_peaks > 1:
                        counter_odd += 2
                    leader_sync_counter.append(counter_odd)
                startvalue = value
        self.sync_counter = leader_sync_counter
        self.is_synced = True
        print("Leader has {} peaks!".format(leader_peaks))
        print("Leader is synced.")    
    
    def loadCalibrationParameters(self, filepath:str):
        """loads  calibrationparameters (A,b) from npz file
        """
       # filepath="Helper/CalibrationParameters/"+filename+".npz"
        import os.path
        if (os.path.isfile(filepath)):
            data=np.load(filepath)
            b=data['b']
            A=data['A']
            return A,b
        else:
          A=np.eye(3)
          b=np.array([0,0,0])
          return A,b
    
    def getCalibrateAccInLSB(self, acc, A, b):
        """ returns Calibrated acc vector
        input is (mx3) acc vector, b = (3x1) sensorbias for acc and A = (3x3) acc Calibration Matrix"""
        accCal=[]
        for a in acc:
            #print(a)
            accCal.append(A@a.T+b.T)
        return np.vstack(accCal)
    
    def calibrate(self, node:str):
        print("Calibrating IMU Data...")
        filepath = ""
        # read in calibration matrix
        if node =="C04":
            filepath = r"C:\git_rep\dataviewer\Calibration\Nodes\C04\C04LAcc070222No1.npz"
        acc_filepath = filepath  #ToDo: extend for gyro, mag
        gyro_filepath = ""
        mag_filepath = ""
        
        # calibrate acc
        A_acc, b_acc = self.loadCalibrationParameters(acc_filepath)
        print("A_acc = \n", A_acc)
        print("b_acc = \n", b_acc)
        acc = self.data.loc[:,'acc_x':'acc_z'].to_numpy()
        acc_cali = self.getCalibrateAccInLSB(acc, A_acc, b_acc)
        self.acc_cali = acc_cali / self.acc_factor * self.gravi
        self.acc_X_cali = self.acc_cali[:, 0]
        self.acc_Y_cali = self.acc_cali[:, 1]
        self.acc_Z_cali = self.acc_cali[:, 2]
        # ToDO: calibrate gyro
        A_gyro, b_gyro = self.loadCalibrationParameters(gyro_filepath)
        gyro = self.data.loc[:,'gyro_x':'gyro_z'].to_numpy()
        gyro_cali = self.getCalibrateAccInLSB(gyro, A_gyro, b_gyro)
        self.gyro_cali = gyro_cali / self.gyro_factor
        self.gyro_X_cali = self.gyro_cali[:, 0]
        self.gyro_Y_cali = self.gyro_cali[:, 1]
        self.gyro_Z_cali = self.gyro_cali[:, 2]
        # ToDo: calibrate mag
        A_mag, b_mag = self.loadCalibrationParameters(mag_filepath)
        mag = self.data.loc[:,'mag_x':'mag_z'].to_numpy()
        mag_cali = self.getCalibrateAccInLSB(mag, A_mag, b_mag)
        self.mag_cali = mag_cali / self.mag_factor
        self.mag_X_cali = self.mag_cali[:, 0]
        self.mag_Y_cali = self.mag_cali[:, 1]
        self.mag_Z_cali = self.mag_cali[:, 2]
        print("Finished Calibration!")


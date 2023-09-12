# -*- coding: utf-8 -*-
"""
Created on 01.03.2023

@author: jfh

Change coordinates for the AvaNode data 

"""


import gps_imu_tools as git
import GPS_Class
from matplotlib import pyplot as plt
import numpy as np


# Load information on AvAnodes data 
gps_L1_C01 = GPS_Class.GPSData()
gps_L2_C01 = GPS_Class.GPSData()
gps_L2_C03 = GPS_Class.GPSData()
gps_L2_C03_2 = GPS_Class.GPSData()
gps_L3_C01 = GPS_Class.GPSData()
gps_L4_C10 = GPS_Class.GPSData()
gps_L5_C07 = GPS_Class.GPSData()
gps_L5_C09 = GPS_Class.GPSData()
gps_L5_C10 = GPS_Class.GPSData()

gps_L1_C01.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L1\ava210315_C01.txt'
gps_L2_C01.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L2\ava210316_C01_GPS.txt'
gps_L2_C03.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L2\ava210316_C03_GPS.txt'
gps_L2_C03_2.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L2\C03Gps002.txt'
gps_L3_C01.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L3\ava220123_C01.txt'
gps_L4_C10.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L4\ava220203_seilbahn_C10_gnss.txt'
gps_L5_C07.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L5\220222_C07_avalanche_GPS.txt'
gps_L5_C09.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L5\220222_C09_avalanche_GPS.txt'
gps_L5_C10.path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\L5\220222_C10_avalanche_GPS.txt'


def processNode(AvaNode,Bez,ax):
    AvaNode.read_data()
    if not AvaNode.vel_bool:
        AvaNode.calc_pos_vel()
    north,east,z = git.gps_to_mercator(AvaNode)
    # AvaNode.plot_vel()
    # AvaNode.plot_vel(plot_pos_vel=True)

    path = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\AvaNode_data\Export' + Bez + r'.txt'
    file = open(path,'w')

    savetotxt(file,north,east,z,AvaNode)
    ax.plot(east,north)


# save coords
def savetotxt(file,north,east,z,AvaNode):
    # print(len(north))
    # file.write('Index,Timestep,North,East,z\n')
    file.write('Time,East,North,z,v\n')

    # Choose Velocity
    print(AvaNode.vel_bool)
    if AvaNode.vel_bool:
        vname = 'v_total'
    else:
        vname = 'v_total_pos'

    for i in range(len(north)):
        # file.write(str(i)+',')
        file.write(str(AvaNode.data.iloc[i][0])+',')
        file.write(str(east[i])+',')
        file.write(str(north[i])+',')
        file.write(str(z[i])+',')
        file.write(str(AvaNode.data_pos[vname][i])+'\n')


fig, ax = plt.subplots()
processNode(gps_L1_C01,r'\L1_C01',ax)
processNode(gps_L2_C01,r'\L2_C01',ax)
processNode(gps_L2_C03,r'\L2_C03',ax)
# processNode(gps_L2_C03_2,r'\L2_C03_preprocessed',ax)          # Hier wurde irgendwas transformiert, für mich erstmal nicht nützlich
processNode(gps_L3_C01,r'\L3_C01',ax)
processNode(gps_L4_C10,r'\L4_C10',ax)
processNode(gps_L5_C07,r'\L5_C07',ax)
processNode(gps_L5_C09,r'\L5_C09',ax)
processNode(gps_L5_C10,r'\L5_C10',ax)

# print(gps_L5_C07.data)
# print(gps_L5_C07.data.iloc[0][0])

# for i in range(10):
#     print(print(gps_L5_C07.data.iloc[i][0]))

# gps_L5_C10.data_pos.keys()[0]




ax.axis('square')
# ax.set_aspect('equal')
ax.set_xlabel('East') 
ax.set_ylabel('North')


plt.show()




# path_c10 = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\AvaNode_data\220222_C10_avalanche_GPS.txt'
# path_c09 = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\AvaNode_data\220222_C09_avalanche_GPS.txt'
# path_c07 = r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Simulationen_Masterarbeit\Vergleich_Nordkette\AvaNode_data\220222_C07_avalanche_GPS.txt'

# gps_c10.path = path_c10
# gps_c09.path = path_c09
# gps_c07.path = path_c07
# gps_c10.read_data()
# gps_c09.read_data()
# gps_c07.read_data()

# n10,e10,z10 = git.gps_to_mercator(gps_c10) 
# n09,e09,z09 = git.gps_to_mercator(gps_c09)
# n07,e07,z07 = git.gps_to_mercator(gps_c07)

# gps_L1_C01.read_data()
# gps_L2_C01.read_data()
# gps_L2_C03.read_data()
# gps_L3_C01.read_data()
# gps_L4_C10.read_data()
# gps_L5_C07.read_data()
# gps_L5_C09.read_data()
# gps_L5_C10.read_data()


# file = open(r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\Export\C07.txt','w')
# savetotxt(file,n07,e07,z07)
# file = open(r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\Export\C09.txt','w')
# savetotxt(file,n09,e09,z09)
# file = open(r'C:\Users\frede\Documents\Programmdateien\AvaFrame\Conda_AvaFrame_GitHub\AvaFrame\avaframe\outPart_JFH\Export\C10.txt','w')
# savetotxt(file,n10,e10,z10)

# ax.plot(e10,n10, color='orange')
# ax.plot(e09,n09, color='green')
# ax.plot(e07,n07, color='brown')
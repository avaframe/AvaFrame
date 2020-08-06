#! /usr/bin/python

import re

def fromString(s):

    #dirty hack to remove newlines in an easy and consistant way
    s = ' '.join(s.splitlines())
    args = s.split(',')

    return fromStringList(args)

def fromStringList(args):

    ad = AimecData()

    for arg in args:
        for attr in ad.__dict__.keys():
            #very dirty hack
            s = re.search(r'.*%s=(.*)'%(attr), arg, re.M|re.I)
            if s:
                if type(ad.__dict__[attr]) is str:
                    ad.__dict__[attr] = s.group(1)
                elif type(ad.__dict__[attr]) is int:
                    ad.__dict__[attr] = int(s.group(1))
                elif type(ad.__dict__[attr]) is bool:
                    ad.__dict__[attr] = (s.group(1).lower in ['true', '1'])
                elif type(ad.__dict__[attr]) is float:
                    ad.__dict__[attr] = float(s.group(1))
                else:
                    raise ValueError('attribute of AimecData from unknown type %s' %(type(ad.__dict__[attr])))

    return ad

class AimecData(object):
    """
    AimecData - parameters for aimec
    """
    #TODO better documentation!!

    def __init__(self):

        self.__explainationDict__ = {}

        self.domainWidth = 500.
        #self.domainLength = 600.
        self.calcPressureLimit = 1.0
#        self.density = 200.
#        self.growthIndex = 3.


        self.pathVelocity = ''
        self.pathFlowHeight = ''
        self.pathPressure = ''
        self.pathEnergy = ''
        self.pathMass = ''
        self.pathNumInfo = ''

        self.pathAvalanchePath = ''
        self.pathDHM = ''
        self.pathDepoArea = ''
#        self.pathEntArea = ''
        self.pathDocDamage = ''
        self.pathAOI = ''
        self.pathDocRadar = ''

        self.pathResult = ''
#        self.pathMapResult = ''


#       explainationDict --> Explanation of fields
        self.__explainationDict__['domainWidth'] = 'width of calculation domain, type integer'
#        self.__explainationDict__['domainLength'] = 'length of calculation domain, type integer'
        self.__explainationDict__['calcPressureLimit'] = 'limit of pressure to calculate run out area, type float'
#        self.__explainationDict__['density'] = 'the density of the mass flow mostly used to convert velocity to pressure'
#        self.__explainationDict__['growthIndex'] = 'the expected (documented) growth index'
        self.__explainationDict__['pathVelocity'] = 'path to the velocity files'
        self.__explainationDict__['pathFlowHeight'] = 'path to the flow height files'
        self.__explainationDict__['pathPressure'] = 'path to the pressure data files'
        self.__explainationDict__['pathEnergy'] = 'path to the enery files'
        self.__explainationDict__['pathMass'] = 'path to something mass related ???'
        self.__explainationDict__['pathNumInfo'] = 'path to nummeric info'
        self.__explainationDict__['pathAvalanchePath'] = 'file with the avalanche path for transformation'
        self.__explainationDict__['pathDHM'] = 'path to the digital simulation elevation model (*.asc with same cellsize like simulation results)'
        self.__explainationDict__['pathDepoArea'] = 'path to the documented run out area'
#        self.__explainationDict__['pathEntArea'] = 'path to the documented entrainment area'
        self.__explainationDict__['pathDocDamage'] = 'path to the documented damage events'
        self.__explainationDict__['pathAOI'] = 'path to area of interest polygon'
        self.__explainationDict__['pathDocRadar'] = 'path to radar input data'
        self.__explainationDict__['pathResult'] = 'path where the aimec results will be written'
#        self.__explainationDict__['pathMapResult'] = 'path where the map results (pdf) will be written, leave empty for no map output'


    def __str__(self):

        string = ''
        for attr in self.__dict__.keys():
            if not '__' in attr:
                if not len(string) == 0: string += ',\n'
                string += attr + '=' + self.__dict__[attr].__str__() + ''

        return string

    def with__(self, arg):

       if len(arg) > 0:
           return True
       else:
           return False


#this is a very basic unittast for the basic IO function of the dataclass
if __name__ == '__main__':
    ad = AimecData()
    ad.domainWidth = 100
    ad.pathPressure = '/home/user/testing/pressure'
    s= ad.__str__()
    ad2 = fromString(s)
    print(ad2.__str__())

import re

def fromString(s):

    #dirty hack to remove newlines in an easy and consistant way
    s = ' '.join(s.splitlines())
    args = s.split(',')

    return fromStringList(args)

def fromStringList(args):

    ad = AvaframeData()

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
                    raise ValueError('attribute of AvaframeData from unknown type %s' %(type(ad.__dict__[attr])))

    return ad

class AvaframeData(object):
    """
    AvaframeData - parameters for avaframe
    """
    #TODO better documentation!!

    def __init__(self):

        self.__explainationDict__ = {}

        self.pathAvalancheName = ''
        self.locCfgAB = 0
        self.locCfg1DFA = 0

#       explainationDict --> Explanation of fields

        self.__explainationDict__['pathAvalancheName'] = 'path where the avalanche folder will be created'
        self.__explainationDict__['locCfgAB'] = '1 if you want to use your custom AB configuration file'
        self.__explainationDict__['locCfg1DFA'] = '1 if you want to use your custom DFA configuration file'

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
    ad = AvaframeData()
    ad.domainWidth = 100
    ad.pathPressure = '/home/user/testing/pressure'
    s= ad.__str__()
    ad2 = fromString(s)
    print(ad2.__str__())

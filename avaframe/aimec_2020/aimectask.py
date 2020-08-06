#! /usr/bin/python

from aimecdata import AimecData

class AimecTask(object):

    def __init__(self):
        self.abort = False

    def cliName(self):
    
        raise NotImplementedError('Calling a pure virtual function')  
        
    def name(self):
        
        raise NotImplementedError('Calling a pure virtual function')
       
    def resultPath(self):
        
        raise NotImplementedError('Calling a pure virtual function')

    def description(self):
        
        raise NotImplementedError('Calling a pure virtual function')
        
    def validateData(self, data):
        
        raise NotImplementedError('Calling a pure virtual function')
        
    def run(self, data, callBack = None):
        
        raise NotImplementedError('Calling a pure virtual function')

#! /usr/bin/python

import aimecdata as aimecdata
from aimecrunner import *


import string
import sys
import re

version_major = 0
version_minor = 0
version_debug = 1
version_build = 0

def versionstring():
    return version_major.__str__() + '.' + version_minor.__str__() + '.' + version_debug.__str__() + '.' + version_build.__str__()

def richtextToTerminal(richtext):
    FG_RED      = '\033[31m'
    FG_GREEN    = '\033[32m'
    FG_YELLOW   = '\033[33m'
    FG_BLUE     = '\033[34m'
    FG_MAGENTA  = '\033[35m'
    FG_CYAN     = '\033[36m'
    FG_DEFAULT  = '\033[39m'

    richtext = richtext.replace('<br/>', '\n')
    richtext = richtext.replace('<b>', FG_YELLOW)
    richtext = richtext.replace('</b>', FG_DEFAULT)
    richtext = richtext.replace('<i>', FG_CYAN)
    richtext = richtext.replace('</i>', FG_DEFAULT)
    richtext = richtext.replace('<em>', FG_MAGENTA)
    richtext = richtext.replace('</em>', FG_DEFAULT)

    return richtext

def main():

    if '--help' in sys.argv:
        usage()
        sys.exit(0)

    if '--v' in sys.argv:
        print(versionstring())
        sys.exit(0)

    #check if there is a file to open
    usefile = False
    for arg in sys.argv:
        s = re.search(r'.*open=(.*)', arg, re.M|re.I)
        if s:
            try:
                f = open(s.group(1), 'r')
            except IOError as e:
                print('I/O error({0}): {1}'.format(e.errno, e.strerror))
                return
            print('using file %s to get parameters' %(s.group(1)))
            aD = aimecdata.fromString(f.read())
            usefile = True
    if not usefile:
        print('reading command line arguments to get parameters')
        aD = aimecdata.fromStringList(sys.argv)

    print('running aimec with dataset:\n' + aD.__str__()+'\n')

    runner = AimecRunner()
    print(runner.tasks)
    #runner.killwhenFinished = True


    tasklist = []
    for task in runner.tasks:
        for arg in sys.argv:
            if arg == task.cliName():
                tasklist.append(task)

    if len(tasklist) == 0:
        print('[Aimec] Error: no task sepcified. Exiting application')
        return
    else:
        print('[Aimec] running the following aimec tasks')
        for task in tasklist:
            print('%s'%(task.name()))
        print('\n')


    #runner.start()
    runner.runall(aD, tasklist)
    #ret = app.exec_()
    del runner
    sys.exit(0)

def usage():
    print('I am aimec %s' %(versionstring()))
    print('usage: %s [options] [parameter=<value>] [task]' %(sys.argv[0]))

    print('\nThe options are:\n')
    print('--help\t\t\t\tprint this message and exit')
    print('--v\t\t\t\tprint version string and exit')
    print('--gui\t\t\t\tlaunch the graphical user interface')
    print('-open=<parameter file>\t\topen the file and read parameters from there. Parameters from the command line get ignored.')

    print('\nThe parameters are:\n')

    ad = AimecData()
    for param in ad.__dict__.keys():
        if not '__' in param:
            print('%s : %s' %(richtextToTerminal('<em>'+param+'</em>'), richtextToTerminal(ad.__explainationDict__[param])))

    print('The parameters are not case sensitive.')
    print('The parameters will be read from command line, from the parameter file (if given) or from the gui.\n')

    print('The following tasks are availible:\n')

    runner = AimecRunner()
    for task in runner.tasks:
        print('%s (%s) \n%s' %(richtextToTerminal('<em>'+task.cliName()+'</em>'), richtextToTerminal('<i>'+task.name()+'</i>'), richtextToTerminal(task.description())))

    print('\n')
    print('\n(c) BFW 2015 - all rights reserved\n')

if __name__ == '__main__':
    main()

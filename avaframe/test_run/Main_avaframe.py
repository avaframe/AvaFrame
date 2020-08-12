import avaframe.test_run.avaframedata as avaframedata
from avaframe.test_run.avaframerunner import *


import string
import sys
import re
import logging

# create logger, set to logging.DEBUG to see all messages
logging.basicConfig(
    level=logging.INFO,
    format="%(module)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("Main_Avaframe.log", "w"),
        logging.StreamHandler()
    ]
)
logmain = logging.getLogger(__name__)

version_major = 0
version_minor = 1


def versionstring():
    return 'Avaframe ' + version_major.__str__() + '.' + version_minor.__str__()


def richtextToTerminal(richtext):
    FG_RED = '\033[31m'
    FG_GREEN = '\033[32m'
    FG_YELLOW = '\033[33m'
    FG_BLUE = '\033[34m'
    FG_MAGENTA = '\033[35m'
    FG_CYAN = '\033[36m'
    FG_DEFAULT = '\033[39m'

    richtext = richtext.replace('<br/>', '\n')
    richtext = richtext.replace('<b>', FG_YELLOW)
    richtext = richtext.replace('</b>', FG_DEFAULT)
    richtext = richtext.replace('<i>', FG_CYAN)
    richtext = richtext.replace('</i>', FG_DEFAULT)
    richtext = richtext.replace('<em>', FG_MAGENTA)
    richtext = richtext.replace('</em>', FG_DEFAULT)
    richtext = richtext.replace('<g>', FG_GREEN)
    richtext = richtext.replace('</g>', FG_DEFAULT)

    return richtext


def main():

    if '--help' in sys.argv:
        usage()
        sys.exit(0)

    if '--v' in sys.argv:
        print(versionstring())
        sys.exit(0)

    # check if there is a file to open
    usefile = False
    for arg in sys.argv:
        s = re.search(r'.*open=(.*)', arg, re.M | re.I)
        if s:
            try:
                f = open(s.group(1), 'r')
            except IOError as e:
                log.error('I/O error({0}): {1}'.format(e.errno, e.strerror))
                return
            log.info('using file %s to get parameters' % (s.group(1)))
            aD = avaframedata.fromString(f.read())
            usefile = True
    if not usefile:
        logmain.info('reading command line arguments to get parameters')
        aD = avaframedata.fromStringList(sys.argv)

    logmain.info('running avaframe with dataset:\n' + aD.__str__()+'\n')

    runner = AvaframeRunner()
    #runner.killwhenFinished = True

    tasklist = []
    for task in runner.tasks:
        for arg in sys.argv:
            if arg == task.cliName():
                tasklist.append(task)


    if not aD.pathAvalancheName:
        logmain.error('No \"pathAvalancheName\" sepcified. Exiting application')
        return
    if len(tasklist) == 0:
        logmain.error('No task sepcified. Exiting application')
        return
    else:
        logmain.info('running the following Avaframe tasks :')
        for task in tasklist:
            logmain.info('\t-%s' % (task.name()))

    # runner.start()
    runner.runall(aD, tasklist)
    #ret = app.exec_()
    del runner
    sys.exit(0)


def usage():
    print('I am avaframe %s' % (versionstring()))
    print('usage: %s [options] [parameter=<value>] [task]' % (sys.argv[0]))

    print('\nThe options are:\n')
    print('--help\t\t\t\tprint this message and exit')
    print('--v\t\t\t\tprint version string and exit')
    print('-open=<parameter file>\t\topen the file and read parameters from there. Parameters from the command line get ignored.')

    print('\nThe parameters are:\n')

    ad = avaframedata.AvaframeData()
    for param in ad.__dict__.keys():
        if not '__' in param:
            print('%s : %s' % (richtextToTerminal('<em>'+param+'</em>'),
                               richtextToTerminal(ad.__explainationDict__[param])))

    print('The parameters are not case sensitive.\n')
    # print('The parameters will be read from command line, from the parameter file (if given) or from the gui.\n')

    print('The following tasks are availible:\n')

    runner = AvaframeRunner()
    for task in runner.tasks:
        print('%s (%s) \n%s\n%s' % (richtextToTerminal('<em>'+task.cliName()+'</em>'), richtextToTerminal('<i>' +
                                                                                                          task.name()+'</i>'), richtextToTerminal(task.description()), richtextToTerminal('<b>'+task.requierments()+'</b>')))

    print('\n')
    print('\n' + richtextToTerminal('<g>'+'(c) BFW 2020 - all rights reserved'+'</g>') + '\n')


if __name__ == '__main__':
    main()

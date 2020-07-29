# !---------written by Felix Oesterle (FSO)------------------
# -DESCRIPTION:
#   Usage: python SHPConv.py import/export FILENAME
#   all modules can be install via pip install
#   shapefile -> pip/pip3 install pyshp
#   example for call
# -----------------------------------------------------------
from __future__ import print_function
import shapefile
import sys
import os
import getopt
import ntpath
import numpy as np


def read_nxyzfile(fname):
    out_array = []

    # read all lines and split at iso's
    with open(fname, 'r') as infile:
        for line in infile:
            if line.startswith("iso"):
                out_array.append(line)
            else:
                out_array[-1] += line

    # FSO--- loop through all entries
    allLines = dict()
    for element in out_array:
        npa = list()
        # FSO--- loop through each line
        for n, entry in enumerate(element.split('\n')):
            if entry.split('=')[0] == 'iso':
                iso = entry.split('=')[1]
            if entry.split('=')[0] == 'unit':
                iso = entry.split('=')[1] + '_' + iso
            if entry.split('=')[0] == 'res':
                iso = entry.split('=')[1] + '_' + iso
            if n > 3:
                floats = [float(x) for x in entry.split()]
                if floats:
                    npa.append(floats[0:2])

        if iso not in allLines:
            allLines[iso] = list()
        allLines[iso].append(npa)

    return(allLines)


def SHP2NXYZ(infile, defname):
    #  The shapefile should have attributes name and d0.
    # FSO--- Input shapefile
    sf = shapefile.Reader(infile)

    # FSO--- Output nxyz file
    outfile = infile+'.nxyz'
    # print('Writing to file: ', outfile)
    nxyz = open(outfile, 'w')

    # FSO--- set defaults for variables
    layername = None
    d0 = None
    rho = None
    sks = None
    iso = None

    # FSO--- get coordinate system
    prjfile = infile.replace('.shp', '.prj')
    if os.path.isfile(prjfile):
        prjf = open(prjfile, 'r')
        sks = prjf.readline()

    # FSO--- Start reading the shapefile
    records = sf.shapeRecords()
    shps = sf.shapes()

    for n, item in enumerate(shps):
        pts = item.points
        zs = [0.0] * len(pts)

        # FSO--- check if records are available and extract
        if records:
            # FSO--- loop through fields
            for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                # FSO--- get entity name
                name = name.lower()
                if (name == 'name'):
                    layername = value
                if (name == 'd0'):
                    d0 = value
                if (name == 'rho'):
                    rho = value
                if (name == 'sks'):
                    sks = value
                if (name == 'iso'):
                    iso = value
            # FSO--- if name is still empty go through file again and take Layer instead
            if ((type(layername) is bytes) or (layername is None)):
                for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                    if (name == 'Layer'):
                        layername = value

        # FSO--- if layer still not defined, use generic
        if layername is None:
            layername = defname+'_'+str(n)

        print('SHPConv: Found layer ', layername)

        # FSO--- Write record header to file
        print('name=', layername, file=nxyz, sep='')

        print('ah=', d0, file=nxyz, sep='')
        print('rho=', rho, file=nxyz, sep='')
        print('sks=', sks, file=nxyz, sep='')

        # FSO--- write each record to file
        print('iso=', iso, file=nxyz, sep='')
        print(len(pts), file=nxyz)
        for (pt, z) in zip(pts, zs):
            print(pt[0], pt[1], z, file=nxyz)


def SHP2Array(infile, defname):
    #  The shapefile should have attributes name and d0.
    # FSO--- Input shapefile
    sf = shapefile.Reader(infile)

    # FSO--- set defaults for variables
    layername = None
    d0 = None
    rho = None
    sks = None
    iso = None

    # FSO--- get coordinate system
    prjfile = infile.replace('.shp', '.prj')
    if os.path.isfile(prjfile):
        prjf = open(prjfile, 'r')
        sks = prjf.readline()

    # FSO--- Start reading the shapefile
    records = sf.shapeRecords()
    shps = sf.shapes()

    SHPdata = {}
    SHPdata['sks'] = sks
    Name = []
    Length = np.empty((0))
    Start = np.empty((0))
    Coordx = np.empty((0))
    Coordy = np.empty((0))
    Coordz = np.empty((0))
    start = 0

    for n, item in enumerate(shps):
        pts = item.points
        zs = [0.0] * len(pts)

        # FSO--- check if records are available and extract
        if records:
            # FSO--- loop through fields
            for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                # FSO--- get entity name
                name = name.lower()
                if (name == 'name'):
                    layername = value
                    Name.append(layername)
                if (name == 'd0'):
                    d0 = value
                if (name == 'rho'):
                    rho = value
                if (name == 'sks'):
                    sks = value
                if (name == 'iso'):
                    iso = value
            # FSO--- if name is still empty go through file again and take Layer instead
            if ((type(layername) is bytes) or (layername is None)):
                for (name, typ, size, deci), value in zip(sf.fields[1:], records[n].record):
                    if (name == 'Layer'):
                        layername = value

        # FSO--- if layer still not defined, use generic
        if layername is None:
            layername = defname+'_'+str(n)

        print('SHPConv: Found layer ', layername)

        Start = np.append(Start, start)
        length = len(pts)
        Length = np.append(Length, length)
        start += length

        for (pt, z) in zip(pts, zs):
            Coordx = np.append(Coordx, pt[0])
            Coordy = np.append(Coordy, pt[1])
            Coordz = np.append(Coordz, z)

    SHPdata['Name'] = Name
    SHPdata['Start'] = Start
    SHPdata['Length'] = Length
    Coord = np.vstack((Coordx, Coordy))
    Coord = np.vstack((Coord, Coordz))
    SHPdata['Coord'] = Coord
    return SHPdata


def NXYZ2SHP(infile, info):
    outfile = infile.replace('.nxyz', '')
    # print('Writing file: ', outfile)

    # FSO--- split info if given, and write prj file
    # if info:
    mdi = dict(it.split("=") for it in info.split(";")[:-2])
    if 'sks' in mdi:
        prjf = open(outfile+'.prj', 'w')
        prjf.write(mdi['sks'])
        prjf.close()

    # FSO--- Input nxyz
    allLines = read_nxyzfile(infile)

    # FSO--- Check whether file is a PSA file (via name)
    # to then output polygon. Default is polyline
    # FSO--- Output shapefile
    psaFlag = False
    if 'psa' in ntpath.basename(infile.lower()):
        print('Seems to be a PSA file')
        psaFlag = True

    # Default to Polyline
    if psaFlag:
        w = shapefile.Writer(outfile, shapeType=5)
    else:
        w = shapefile.Writer(outfile, shapeType=3)

    w.autoBalance = 1
    w.field('Layer', 'C', '40')
    w.field('Info', 'C', '40')

    # Make sure all lines are closed
    for key in allLines:
        print(key)
        for line in allLines[key]:
            if not line[0] == line[-1]:
                line.append(line[0])

    for key in allLines.keys():
        # write line parts
        if psaFlag:
            w.poly(allLines[key])
        else:
            w.line(allLines[key])

        # write the fields
        w.record(key, 'From nxyz')

    # w.save(outfile)
    w.close()

    os.remove(infile)


if __name__ == "__main__":

    defname = 'SHP'
    info = 'info'

    # FSO--- get Inputfile and action from commandline arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:a:n:s:',
                                   ['infile=', 'action=', 'name=', 'sks='])
    except getopt.GetoptError:
        print('SHPConv.py -i <inputfile> -a <import/export>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('SHPConv.py -i <inputfile> -a <import/export>')
            sys.exit()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-a", "--action"):
            action = arg
        elif opt in ("-n", "--name"):
            defname = arg
        elif opt in ("-s", "--sks"):
            info = arg

    print('SHPConv: ', action, ' file: ', infile)

    # import converts shp -> nxyz
    if action == 'import':
        SHP2NXYZ(infile, defname)
    # import converts nxyz -> shp
    elif action == 'export':
        NXYZ2SHP(infile, info)
    else:
        print('SHPConv: Need import or export as -a cmdline argument')

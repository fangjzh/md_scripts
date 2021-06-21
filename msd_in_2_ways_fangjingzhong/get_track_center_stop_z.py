#!/usr/bin/python3
import numpy as np
import math
import sys


stop_flag = 0
stop_val = 80
timestep_a = 0
## the file is not large, so we use readlines, fdata is a list not a iter
## if the file is large, readline is recommended
readfile = sys.argv[1]
timestep = 0
fw = open('atomic.dat','w')
with open(readfile,'r') as fh:
    line = fh.readline()
    while line:
        line = fh.readline()
        timestep_b = timestep
        timestep = int(line)
        line = fh.readline()
        line = fh.readline()
        num_atom = int(line)
        line = fh.readline()
        line = fh.readline()
        [xl, xh] = map(float,line.split()[0:2])
        line = fh.readline()
        [yl, yh] = map(float,line.split()[0:2])
        line = fh.readline()
        [zl, zh] = map(float,line.split()[0:2])
        line = fh.readline()
        lenx = (xh - xl)/2
        leny = (yh - yl)/2
        lenz = (zh - zl)/2
        sumx = 0.0
        sumy = 0.0
        sumz = 0.0
        for i in range(num_atom):
            line = fh.readline()
            [x, y, z] = map(float,line.split()[0:3])
            if (i == 0):
                x0 = x
                y0 = y
                z0 = z
            if (x - x0 > lenx):
                x = x - 2*lenx
            if (x - x0 < -lenx):
                x = x + 2*lenx
            if (y - y0 > leny):
                y = y - 2*leny
            if (y - y0 < -leny):
                y = y + 2*leny
            if (z - z0 > lenz):
                z = z - 2*lenz
            if (z - z0 < -lenz):
                z = z + 2*lenz
            sumx = sumx + x
            sumy = sumy + y
            sumz = sumz + z
        sumx = sumx / num_atom
        sumy = sumy / num_atom
        sumz = sumz / num_atom
        if (sumx < xl):
            sumx = sumx + 2*lenx
        if (sumx > xh):
            sumx = sumx - 2*lenx
        if (sumy < yl):
            sumy = sumy + 2*leny
        if (sumy > yh):
            sumy = sumy - 2*leny
        if (sumz < zl):
            sumz = sumz + 2*lenz
        if (sumz > zh):
            sumz = sumz - 2*lenz
        if (sumz < -stop_val or sumz > stop_val):
            if (stop_flag == 0):
                timestep_a = timestep
                print("defects is out of the fence at %8d\n" %(timestep))
                fw.write("STEP:   %8d\n" % (timestep))
                fw.write("%f %f %f\n" % (sumx,sumy,sumz))
                stop_flag = 1
        if (stop_flag == 0):
            fw.write("STEP:   %8d\n" % (timestep))
            fw.write("%f %f %f\n" % (sumx,sumy,sumz))
        line = fh.readline()

if (timestep_a == 0):
    timestep_a = timestep

fw.close()

fw2 = open("prog_in.txt",'w')
fw2.write("Control parametre file\nAtom: He  Tempreture: 100K\nBox size in angstrom\n")
fw2.write("%f %f %f\n" % (lenx*2,leny*2,lenz*2))
fw2.write("Number of atoms:\n")
fw2.write("%8d\n" % (1))
fw2.write("Total number of data sets:\n")
fw2.write("%8d\n" % (int(timestep_a/(timestep-timestep_b))))
fw2.write("Interval time in fs (between two sets of data)\n")
fw2.write("%8d\n" % (int(timestep-timestep_b)))
fw2.close()

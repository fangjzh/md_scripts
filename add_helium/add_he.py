#!/bin/env python3

import numpy as np
import math
import sys
import random

### add atoms/total*1000000
appm = 685
bnx = 30
bny = 30
bnz = 160

def gen_rand_atoms(num_tot,num_add):
    num_flag = [0]*num_tot
    while (num_add > 0):
        iseed = random.randint(1,num_tot)
        if (num_flag[iseed - 1] != 1):
            num_flag[iseed - 1] = 1
            num_add = num_add - 1
    return num_flag


## the file is not large, so we use readlines, fdata is a list not a iter
## if the file is large, readline is recommended
#readfile = sys.argv[1]
readfile = 'G:\\Caculation_Data\\temperature_gradient\\W-He\\build_box\\rd.dat'

fw = open('G:\\Caculation_Data\\temperature_gradient\\W-He\\build_box\\rd_r.dat','w')

with open(readfile,'r') as fh:
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    num_atom = int(line.split()[0])
    num_add = int(num_atom/1000000.0*appm)
    num_flag = gen_rand_atoms(num_atom,num_add)
    fw.write(str(num_atom+num_add) + ' atoms\n')
    line = fh.readline()
    fw.write(line)
    num_type = int(line.split()[0])
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    fw.write(line)
    [xl, xh] = map(float,line.split()[0:2])
    line = fh.readline()
    fw.write(line)
    [yl, yh] = map(float,line.split()[0:2])
    line = fh.readline()
    fw.write(line)
    [zl, zh] = map(float,line.split()[0:2])
    lenx = xh - xl
    a0 = lenx/bnx
    leny = yh - yl
    lenz = zh - zl
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    fw.write(line)
    if(line.split()[0] == "Masses"):
        line = fh.readline()
        fw.write(line)
        mass = [0.0]*num_type
        for i in range(num_type):
            line = fh.readline()
            fw.write(line)
            mass[i] = float(line.split()[1])
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    fw.write(line)
    if(line.split()[0] == "Atoms"):
        line = fh.readline()
        fw.write(line)
        iadd = 0  # current number of added atoms
        for i in range(0,num_atom):
            line = fh.readline()
            fw.write(line)
            atom_id = int(line.split()[0])
            atom_type =  line.split()[1]
            if(num_flag[atom_id - 1] == 1):
                iadd = iadd + 1
                [x,y,z] = map(float,line.split()[2:5])
                x = x + a0*0.25
                if(x > xh):
                    x = x - xh + xl
                y = y + a0*0.5
                if(y > yh):
                    y = y - yh + yl
                fw.write(str(num_atom + iadd)+" 2 "+str(x)+" "+str(y)+" "+str(z)+" 0 0 0\n") 
    line = fh.readline()
    fw.write(line)
    line = fh.readline()
    fw.write(line)
    if(line.split()[0] == "Velocities"):
        line = fh.readline()
        fw.write(line)
        iadd = 0  # current number of added atoms
        for i in range(0,num_atom):
            line = fh.readline()
            fw.write(line)
            atom_id = int(line.split()[0])
            if(num_flag[atom_id - 1] == 1):
                [vx,vy,vz] = map(float,line.split()[1:4])
                fw.write(str(atom_id) +" "+ str(vx) + " " +str(vy) +" "+ str(vz) +"\n")    

fw.close()




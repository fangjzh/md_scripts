#!/usr/bin/env python3
## write by Fangjzh@foxmail.com ####
## 2021-06-20 first edt. ###

import numpy as np

###range###
r_l = -3
r_h = 3


print("Solve the orthorhombic crystal orientation, \
the integer solution range is [%d,%d],\
    original lattice vector is [1 0 0] [0 1 0] [0 0 1]." % (r_l,r_h))
r_h = r_h + 1

line=input("please input an intager vector as x oritation: ")
a = line.split()
a = np.array([ int(x) for x in a ])
if ( len(a) != 3 ):
    print("error input!!! exit!")
    exit()
b = np.empty(shape=[0,3], dtype=int)
for i in range(r_l,r_h):
    for j in range(r_l,r_h):
        for k in range(r_l,r_h):
            mod2_r = np.linalg.norm(np.array([i,j,k])%2)
            mod3_r = np.linalg.norm(np.array([i,j,k])%3)
            if (mod2_r != 0 and mod3_r != 0):
                if ( a.dot(np.array([i,j,k]) ) == 0 ):
                    b = np.append(b,[[i,j,k]],axis=0)

if ( len(b) == 0 ):
    print("no result in [-3,3]!!! exit!")
    exit()
num_r = 0

for b1 in b:
    for i in range(r_l,r_h):
        for j in range(r_l,r_h):
            for k in range(r_l,r_h):
                mod2_r = np.linalg.norm(np.array([i,j,k])%2)
                mod3_r = np.linalg.norm(np.array([i,j,k])%3)
                if (mod2_r != 0 and mod3_r != 0):
                    cross_r = np.cross(a,b1)
                    ax = np.array([i,j,k])/np.linalg.norm(np.array([i,j,k]))
                    bx = cross_r/np.linalg.norm(cross_r)
                    cross_r1 = np.linalg.norm(ax + bx)  ## If the two vectors are in opposite directions 
                                                        ## cross_r1 is 0, must be very small
                    cross_r = np.cross(cross_r,np.array([i,j,k]))
                    cross_r2 = np.linalg.norm(cross_r)
                    if ( cross_r2 == 0  and cross_r1 > 0.1):
                        num_r = 1
                        print(a,b1,[i,j,k])

if (num_r == 0):
    print("no result!!")


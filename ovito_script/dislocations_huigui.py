#!/usr/bin/env python
# -- coding: utf-8 --
#### 对于位错环的惯性面计算似乎不准确，虽然回归好像问题不大
import sys
import math
import numpy as np

from sklearn import linear_model

in_file = sys.argv[1]
out_file = "result.txt"
o_file = open(out_file,'w')

atom_file = open(in_file)

xx = [0.0]*60000
yy = [0.0]*60000
zz = [0.0]*60000


i = 0
a_num = 0

for line in atom_file:
  i = i + 1
  atom = line.split()
  if i == 5:
    a_num = int(atom[1])
  elif i > 5:
    if i < a_num + 6 :
      xx[i-6] = float(atom[0])
      yy[i-6] = float(atom[1])
      zz[i-6] = float(atom[2])

atom_file.close()

# 构建成特征、值的形式
X,Z = np.column_stack((np.array(xx[0:a_num]),np.array(yy[0:a_num]))), np.array(zz[0:a_num])

# 建立线性回归模型
regr = linear_model.LinearRegression()

# 拟合
regr.fit(X, Z)


# 不难得到平面的系数、截距
a, b = regr.coef_, regr.intercept_

a = np.append(a,-1.0)
la = np.sqrt(a.dot(a))
a = a/la
print('normal vector of the surface:')
print(a)



o_file.write('normal vector of the surface:'+'\n')
o_file.write(str(a)+'\n')



###########
nr = np.array([1.0,1.0,1.0])


lnr = np.sqrt(nr.dot(nr))

cos_angle = a.dot(nr)/(la*lnr)
angle = np.arccos(cos_angle)
angle2 = angle*180/np.pi

if(angle2 > 90):
  nr = 0 - nr
  angle2 = 180 - angle2
print(nr)
print("the angle is: %f" % (angle2) )
o_file.write(str(nr) + '\n')
o_file.write(str(angle2) + '\n')
############

###########
nr = np.array([1.0,-1.0,1.0])


lnr = np.sqrt(nr.dot(nr))

cos_angle = a.dot(nr)/(la*lnr)
angle = np.arccos(cos_angle)
angle2 = angle*180/np.pi
if(angle2 > 90):
  nr = 0 - nr
  angle2 = 180 - angle2
print(nr)
print("the angle is: %f" % (angle2) )
o_file.write(str(nr) + '\n')
o_file.write(str(angle2) + '\n')
############

###########
nr = np.array([1.0,-1.0,-1.0])


lnr = np.sqrt(nr.dot(nr))

cos_angle = a.dot(nr)/(la*lnr)
angle = np.arccos(cos_angle)
angle2 = angle*180/np.pi
if(angle2 > 90):
  nr = 0 - nr
  angle2 = 180 - angle2
print(nr)
print("the angle is: %f" % (angle2) )
o_file.write(str(nr) + '\n')
o_file.write(str(angle2) + '\n')
############

###########
nr = np.array([-1.0,-1.0,1.0])


lnr = np.sqrt(nr.dot(nr))

cos_angle = a.dot(nr)/(la*lnr)
angle = np.arccos(cos_angle)
angle2 = angle*180/np.pi
if(angle2 > 90):
  nr = 0 - nr
  angle2 = 180 - angle2
print(nr)
print("the angle is: %f" % (angle2) )
o_file.write(str(nr) + '\n')
o_file.write(str(angle2) + '\n')
############


o_file.close()




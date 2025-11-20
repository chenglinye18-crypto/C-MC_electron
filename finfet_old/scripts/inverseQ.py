#!/usr/bin/python 
#coding=utf-8

import sys, re, os

#read lgrid.txt to get mesh

def read_grid():

  file = open('../lgrid.txt')

  len = int(file.readline())
  counter = 0
  x = []
  while (counter < len) :
    x.append(float(file.readline()) * 1000)
    counter += 1

  len = int(file.readline())
  counter = 0
  y = []
  while (counter < len) :
    y.append(float(file.readline()) * 1000)
    counter += 1

  len = int(file.readline())
  counter = 0
  z = []
  while (counter < len) :
    z.append(float(file.readline()) * 1000)
    counter += 1

  file.close()
  
  return [x,y,z]

#find index according to value

def index(l, val):
  counter = 0
  while ((counter < len(l)) and (abs(l[counter] - val) > 1e-6)) :
    counter += 1
  if (counter >= len(l)):
    raise ValueError
  return counter

#read volume file 

def read_volume(x, y, z):
  file = open('../data/pvolume')

  done = 0
  v = [[[0 for raw in range(len(z))] for raw in range(len(y))] for raw in range(len(x))]
  while not done: 
    line = file.readline()
    if line == "":
      done = 1
    else :
      line_list = re.split('\s+', line)

      i = int(line_list[0])
      j = int(line_list[1])
      k = int(line_list[2])

      v[i][j][k] = float(line_list[3])

  file.close();
  return v
    
def calculate_charge(x, y, z, xr, yr, zr, v, name):

  file = open(name)
  done = 0
  q = 0
  vol = 0

  while not done: 
    line = file.readline()
    if line == "":
      done = 1
    else :
      line_list = re.split('\s+', line)

      i = int(line_list[0])
      j = int(line_list[1])
      k = int(line_list[2])

      if ((i >= xr[0]) and (i <= xr[1]) and (j >= yr[0]) and (j <= yr[1]) and (k >= zr[0]) and (k <= zr[1])):
	dens = float(line_list[7])
	q += dens * v[i][j][k]
        vol += v[i][j][k]

  file.close()

  return q / vol
 
#main funciton

if (len(sys.argv) != 7) :
  print 'error: the ranges of xyz direction are needed'
  sys.exit()
else :
  print '\nthe range of channel is [%s,%s]' % (sys.argv[1], sys.argv[2])
  print 'range of x: [%s,%s]' % (sys.argv[3], sys.argv[4])
  print 'range of z: [%s,%s]' % (sys.argv[5], sys.argv[6])

[x,y,z] = read_grid()

yr = [index(y, float(sys.argv[1])), index(y, float(sys.argv[2]))]
  
xr = [index(x, float(sys.argv[3])), index(x, float(sys.argv[4]))]

zr = [index(z, float(sys.argv[5])), index(z, float(sys.argv[6]))]

v = read_volume(x, y, z)

files = os.listdir('../data')

counter = 0
steps = []
for fname in files:
  if re.match('Electron', fname):
    counter += 1
    steps += [int(fname[len('Electron'):])]

steps.sort()

Q = 0
print
for step in steps[1:]:
  name = '../data/Electron' + str(step)
  Qstep = calculate_charge(x,y,z, xr, yr, zr, v, name)
  Q += Qstep
  print 'charge of step %s : %E' % (step, Qstep)

print 'the average of inverse charge in channel is : %E' % (Q / (len(steps) - 1))


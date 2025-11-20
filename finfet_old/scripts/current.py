#!/usr/bin/python
#coding=utf-8
import re
#打开保存电流的文件
file = open('../data/current')
#读文件
allLines = file.readlines()
file.close()

steps = 0
cur = 0
icont = 0
items = 1

for eachLine in allLines:
  if re.match('step', eachLine):
     steps = steps + 1
  if re.match('icont', eachLine):
    list = re.split('\s+', eachLine)
    icont = int (list[2])
#不用前 items 次统计的数据
  if ((icont <= 1) and (steps > items) and (re.match('Current', eachLine))):
    list = re.split('\s+', eachLine)
    cur = cur + abs(float(list[2]))

print "total steps = %d, last %d items used, and the current is %E" % (steps, steps - items, cur / (2 * (steps - items)))

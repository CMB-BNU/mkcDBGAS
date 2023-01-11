#!/usr/local/bin/python
# -*- coding: utf-8 -*-
#coding=utf-8
from sys import argv
unique=open(argv[1]+".unique","w")
ioe = open(argv[1],"r")
done = {}
while 1:
    lines = ioe.readlines(100000)
    if not lines:
        break
    for line in lines:
        line = line.strip('\n')
        hang = line.split('\t')
        queryName = hang[0]
        subjectName = hang[1]

        bothname=""
        if queryName < subjectName:
            bothname = queryName+","+subjectName
        elif queryName > subjectName:
            bothname = subjectName+","+queryName
        #print(bothname)
        if queryName != subjectName and bothname not in done:#将不是与自己相比较的去除掉，去除hsp为空的情况
            done[bothname] = 1
            print(bothname,file=unique)
t_max = int(argv[2])
print(len(done))
cut_order = int(len(done)/t_max)+1
num = 0
t = 1
uuu=open(argv[1]+".split"+str(t),"a")
for name in done:
	num = num + 1
	cut = cut_order * t
	if num < cut:
		print(name,file=uuu)
	else:
		print(name,file=uuu)
		t = t + 1
		uuu.close()
		uuu=open(argv[1]+".split"+str(t),"a")


#print(datetime.datetime.now())

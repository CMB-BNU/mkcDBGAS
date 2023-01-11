#!/usr/local/bin/python
# -*- coding: utf-8 -*-
#coding=utf-8

import re
import build
from Bio import SeqIO
from sys import argv
import datetime

#print(datetime.datetime.now())
seqs={}
#for seq_record in SeqIO.parse("Araport11_cdna_onlyname.fasta", "fasta"):
for seq_record in SeqIO.parse(argv[2], "fasta"):
	seqs[seq_record.id]=seq_record.seq._data

def secondbubble_old(bubble,queryName,subjectName):
	bubble2={}
	bu_number=len(bubble)+1
	for x in bubble:
		if bubble[x]["length2"] > 5 and  bubble[x]["length1"] > 5:
			bbbb=build.secondconstruteDG(bu_number,seqs[queryName][int(bubble[x]["start_pos1"]):int(bubble[x]["end_pos1"])],seqs[subjectName][int(bubble[x]["start_pos2"]):int(bubble[x]["end_pos2"])])
			bu_number=bu_number+len(bbbb)
			bubble2.update(bbbb)
	return bubble2

def secondbubble(bubble,seq1,seq2):
	bubble2={}
	bubble2.update(bubble)
	for x in bubble:
		if bubble[x]["length2"] > 5 and  bubble[x]["length1"] > 5:
			bbbb=build.secondconstruteDG(bubble[x]["start_pos1"],bubble[x]["start_pos2"],x,seq1[int(bubble[x]["start_pos1"])-1:int(bubble[x]["end_pos1"])],seq2[int(bubble[x]["start_pos2"])-1:int(bubble[x]["end_pos2"])])
			bubble2.update(bbbb)#bbbb是bubble1.1 bubble1.2等key
	return bubble2
def secondbubble(bubble,seq1,seq2):
	bubble2={}
	bubble2.update(bubble)
	for x in bubble:
		if bubble[x]["length2"] > 5 and  bubble[x]["length1"] > 5:
			bbbb=build.secondconstruteDG(bubble[x]["start_pos1"],bubble[x]["start_pos2"],x,seq1[int(bubble[x]["start_pos1"])-1:int(bubble[x]["end_pos1"])],seq2[int(bubble[x]["start_pos2"])-1:int(bubble[x]["end_pos2"])])
			bubble2.update(bbbb)#bbbb是bubble1.1 bubble1.2等key
	return bubble2

af=open(argv[1]+".AF","w")
al=open(argv[1]+".AL","w")
mx=open(argv[1]+".MX","w")
snp=open(argv[1]+".snp","w")
bb=open(argv[1]+".bubble","w")

ioe = open(argv[1],"r")
done=open(argv[1]+".done","w")

while 1:
	lines = ioe.readlines(100000)
	if not lines:
		break
	for line in lines:
		line=line.strip('\n')
		hang = line.split(',')
		queryName = hang[0]
		subjectName = hang[1]
		bothname = queryName+","+subjectName
		print(bothname,file=done)
		if queryName != subjectName:#将不是与自己相比较的去除掉
			bubble=build.construteDG(seqs[queryName],seqs[subjectName])#queryName is seq1
			bubble=secondbubble(bubble,seqs[queryName],seqs[subjectName])
			flag = "0"
			snp_num=0
			outline = ''
			outline_af = ''
			outline_al = ''
			outline_mx = ''
			bridge = 0#coverage
			lastendpos = 0
			if len(bubble) == 1 and bubble["bubble1"]["length2"] > 10 and bubble["bubble1"]["length1"] > 10:
				#print(bubble)
				if bubble["bubble1"]["length1"] / len(seqs[queryName]) <= 0.3 or bubble["bubble1"]["length2"] / len(seqs[subjectName]) <= 0.3:
					if "Q" in bubble["bubble1"]["start_node"]:
						flag = "af"
						outline_af = bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "\n"
					elif "Q" in bubble["bubble1"]["end_node"]:
						flag = "al"
						outline_al = bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "\n"
					else:
						flag = "mx"
						outline_mx = bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "\n"
					bridge = len(seqs[queryName]) - bubble["bubble1"]["length1"]
			else:
				for x in bubble:
					bridge = bridge + int(bubble[x]["start_pos1"]) - lastendpos
					lastendpos = int(bubble[x]["end_pos1"])
					if bubble[x]["length1"] == 0 and bubble[x]["length2"] > 2 and "Q" not in bubble[x]["start_node"] and "Q" not in bubble[x]["end_node"]:#length2为0，说明第一个序列有，所以是query的
						flag = "bubble"+"_"+str(bubble[x]["length1"])+"_"+str(bubble[x]["length2"])
						outline = outline + bothname + "," + seqs[subjectName][max(0,int(bubble[x]["start_pos2"]) - 50):int(bubble[x]["start_pos2"])] + "," + seqs[subjectName][int(bubble[x]["start_pos2"]):int(bubble[x]["end_pos2"])] + "," + seqs[subjectName][int(bubble[x]["end_pos2"]):int(bubble[x]["end_pos2"])+50] + ',' + str(bubble[x]["start_pos1"]) + ',' + str(bubble[x]["end_pos1"]) + ',' + str(bubble[x]["start_pos2"]) + ',' + str(bubble[x]["end_pos2"]) +','+ str(len(seqs[queryName])) +','+ str(len(seqs[subjectName])) + "\n"
					elif bubble[x]["length2"] == 0 and bubble[x]["length1"] > 2 and "Q" not in bubble[x]["start_node"] and "Q" not in bubble[x]["end_node"]:#length2为0，说明第一个序列有，所以是query的
						flag = "bubble"+"_"+str(bubble[x]["length1"])+"_"+str(bubble[x]["length2"])
						outline = outline + bothname + "," + seqs[queryName][max(0,int(bubble[x]["start_pos1"]) - 50):int(bubble[x]["start_pos1"])] + "," + seqs[queryName][int(bubble[x]["start_pos1"]):int(bubble[x]["end_pos1"])] + "," + seqs[queryName][int(bubble[x]["end_pos1"]):int(bubble[x]["end_pos1"])+50] + ',' + str(bubble[x]["start_pos1"]) + ',' + str(bubble[x]["end_pos1"]) + ',' + str(bubble[x]["start_pos2"]) + ',' + str(bubble[x]["end_pos2"]) +','+ str(len(seqs[queryName])) +','+ str(len(seqs[subjectName])) + "\n"
					elif bubble[x]["length2"] == 1 and bubble[x]["length1"] == 1:
						snp_num=snp_num+1
					elif bubble[x]["length2"] > 10 and bubble[x]["length1"] > 10:
						if "Q" in bubble["bubble1"]["start_node"]:
							flag = "af"
							outline_af = outline_af + bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "," + "\n"
						elif "Q" in bubble["bubble1"]["end_node"]:
							flag = "al"
							outline_al = outline_al + bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "," + "\n"
						else:
							flag = "mx"
							outline_mx = outline_mx + bothname + "," + str(bubble["bubble1"]["start_pos1"]) + "," + str(bubble["bubble1"]["end_pos1"]) + "," + str(bubble["bubble1"]["start_pos2"]) + "," + str(bubble["bubble1"]["end_pos2"]) + "," + "\n"
			bridge = bridge + len(seqs[queryName])-lastendpos
			coverage = bridge / min(len(seqs[queryName]), len(seqs[subjectName]))
			if flag != "0" and snp_num > 2 and coverage > 0.6:
				print(snp_num,outline,end='',file=snp)
			if flag != "0" and snp_num < 3 and coverage > 0.6:
				if ",," not in outline:
					print(outline,end='')
				print(outline_af,end='',file=af)
				print(outline_al,end='',file=al)
				print(outline_mx,end='',file=mx)
				for x in bubble:
					print(queryName+","+subjectName,end=",",file=bb)
					print(x,bubble[x],file=bb)

				print(queryName+","+subjectName+","+flag+","+str(snp_num),file=bb)
#print(datetime.datetime.now())

#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#
import re
from sys import argv

def features(seq):
	length = len(seq)
	structure(seq,length)
	fraquency(seq,length)
	distribution(seq,length)
	splitsite(seq)
	
def structure(seq,length):
	#三的倍数、GC含量
	t = (1 if(length%3 == 0) else 0)
	count_GC = len(re.findall('[GCgc]', seq))
	print(t,count_GC,end=",",sep=",")

	#终止密码子
	for s in stop:
		print(seq.count(s),end=",",sep=",")

#单核苷酸，二连核苷酸，三联核苷酸频率，4+16+64=84个特征
def fraquency(seq,length):
	for f in nucl:
		print(round(seq.count(f)/length,2),end=",",sep=",")
		for s in nucl:
			if length > 1:
				print(round(seq.count(f+s)/(length-1),2),end=",",sep=",")
			else:
				print(0,end=",",sep=",")
			for t in nucl:
				if length > 2:
					print(round(seq.count(f+s+t)/(length-2),2),end=",",sep=",")
				else:
					print(0,end=",",sep=",")
				

#四个核苷酸在序列中的分布，第一个，25%，50%，75%，最后一个的位置。20个特征
def distribution(seq,length):
	for x in nucl:
		index_list = [i.start()+1 for i in re.finditer(x, seq)]
		if index_list:
			print(index_list[0]/length,index_list[int(len(index_list)/4)-1]/length,index_list[int(len(index_list)/2)-1]/length,index_list[int(len(index_list)/4*3)-1]/length,index_list[int(len(index_list)/4*4)-1]/length,end=",",sep=",")
		else:
			print(0,0,0,0,0,end=",",sep=",")

	
def splitsite(seq):
	splitsite=("GT","GC","AT","AG","AC")
	for m in splitsite:
		doner = seq[-5:]
		acceptor = seq[:5]
		print(doner.count(m),end=",",sep=",")
		print(acceptor.count(m),end=",",sep=",")

def Dmotif(seq):
	dmotif=("GTAAC","AAGTGT","GTTTGT","ATTAACA","TGAAG","TAACC","TTGAAAT","AATTG","CTGCT","TTTATG","TGATAAA","ATGTTT","TTTCCAA","AAGTC","AAAGA","GTACGT","GTTAAA","GAGCTG","TCATTTT","TGCATG","TTCTT","TTTATC","ACATTT","TGCCAGC","ATAATT","GTAGG","GTATCCT","CATTTG","ACTAAC","TTTCAG","AATTGA","TTAGCA","CAAAT","TAATG","TTTTGAT","AGAAAT","TTTCTA","TATTTC","TAACT","TGAGG","TAAAAT","TTTATA","TATCCT","GTTAGT","TTTACAG","TATTTG","GTACTGT","TTAAG","CATAAA","GCATG","TGATTA","TTTTAAA","CTGACT","ACTAAT","GAGTA","TCTTAA","TTGGTT","ATATTT","AGAGCCA","TCTTT","AGTTTT","GTATTT","TCAGA","TAAGT","AAGCA","TTCACAG","AGTAA","TCTGG","AGCTTT","TGATTTG","TTTTGC","TAGAAA","GTGAG","TTCTGT","GTAAG","ATGAAA","AGAAAA","TGAGC","TGGCTT","TTAATCT","AATTAT","TGGAAAT","CCACAG","AAATGA","GCAAGT","GTAAAA","GTCTG","AAATGT","TGCAT","GAGAAA","TTAGA","TTTATAA","GTTTT","GCTTGGC","TAAGC","GTATG","AAATT","GTAAT","TTCTCT","TGAGAA","TTAGTT","TAAGG","TGTTTAA","GTCAGT","AGAATT","TAAATG","AATTCA","TCCTTT","TAAGA","AAATCA","TAATTTG","GAAATA","TGGTTT","TGTTAA","TGTCT","GTTGGT","TGAATT","AATTTA","TATGT","GCATTT","AAGTA","GCTTCT","TTCTAA","GTTTCT","AGATTT","GAAAAT","TGCTAA","AAGCT","CTTTGCT","TCTGA","TTTCTC","TTTATTC","TTTGCC","TGAAAG","TGTTCT","CTTTT","TTTTCTG","TGAGT","TTGCAG","TAATA","AGTAT","ATTCT","TGCCTTT","ATCAAA","GAGTG","TAGGT","CTTTA","TTTAG","TGATTTT","TTTCAT","CTTTCA","AAGAT","TGCTT","GTGGGT","GTAAAG","CTGAA","TCTGC","CTAAA")
	for m in dmotif:
		if m in seq:
			print(1,end=",",sep=",")
		else:
			print(0,end=",",sep=",")
def Umotif(seq):
	umotif=("GTTTGT","TCTCC","GATTTT","TTTTTC","TAACC","TTGAAAT","AAGCCA","AATTG","CTGCT","TTTATG","TTCACA","TGATAA","ATGTTT","TCCAG","TTTCCAA","TTATTTC","TGTGTT","TCTTG","TTGTAA","CTTGAC","TTAAAAC","CTAAC","AAAGCT","TCTTC","TGCATG","TTCTT","ACATTT","TTTATC","ATTTTCT","GCTGACC","ATAATT","AACAG","CATTTG","TTTCAG","TTAGCA","TTGCCT","CAAAT","TTTAAC","TAATG","AGAAAT","TTTTGAT","TAACT","TAAAAT","AATTACA","TTCAAAA","TTTATA","CTTGTC","TTTACAG","TGGATTT","TTGCATT","TTAAG","TTGGT","CATAAA","TGATTA","GCTTTGC","ATTAG","TTTTAAA","ACTAAT","CTGACT","ATATTT","CATTTA","AAATCT","TTTTGGC","TCTTT","TTATTGA","TCAGA","TAAGT","CTCTG","TCTGG","AATTC","TTTTCC","TTTTGC","GTGAG","TCCATTT","AATTTT","CTTGATT","GTAAG","ATGAAA","AGAAAA","TGGCTT","CTCAG","TGGAAAT","AATTAT","AATAAT","TCCTAG","CCACAG","TCATTTC","AAAGCA","AAATGA","TTTATAG","ATTAAAT","CCTGCAG","TTACAG","AAATGT","TGCAT","CTTCT","TTAGAA","TGTTTC","TTTAC","GTTTT","CTTCCA","TTCTAG","AAATT","TTAAAC","TGAGAA","TTTGTAG","GTCAGT","TAAGA","AAATCA","TGTTGA","CTTGC","TAATTTG","CCTCT","TGGTTT","TGATTTC","TGTTAA","TGTGTC","TCTCT","TTAACA","TTTGGT","TGTCT","TTCCTT","TGAATT","AATTTA","GTTTCT","TGCTAA","AGATTT","GAAAAT","GTTTAAT","TTTGACT","TCTGA","TCTGTT","TTTATTC","TGAAAG","TGTTCT","CTTTT","ATTTGT","TGAGT","CCCCAG","TTGCAG","CTGAT","TAATA","TCTTA","ATTCT","ATCAAA","CTTTA","TTTAG","TTGCTG","TTTCAT","CTTTCA","TTCTC","TGCTT","TCTGC","CTGAA","GTAGGT","CTAAA")
	for m in umotif:
		if m in seq:
			print(1,end=",",sep=",")
		else:
			print(0,end=",",sep=",")
#定义变量
nucl=("A","T","C","G")
stop=("TAA","TAG","TGA")

ioe = open(argv[1],"r")
event = open(argv[2],"w")

while 1:
    lines = ioe.readlines(100000)
    if not lines:
        break
    for line in lines:
        line=line.strip('\n')
        line=line.replace("a","A").replace("g","G").replace("c","C").replace("t","T")#都换成大写
        seqs = line.split(",")
        up = seqs[2]
        a3 = seqs[3]
        down=seqs[4]
        bothname = seqs[0]+","+seqs[1]+","+seqs[5]+","+seqs[6]+","+seqs[7]+","+seqs[8]
        print(bothname,file=event)
        print(len(a3),end=",",sep=",")
        #1
        Dmotif(a3+down)#是否有motif
        Umotif(up+a3)#是否有motif
        #2
        features(up+a3+down)
        #3
        features(a3)
        #4
        features(up)
        features(down)
        #5
        features(up[-30:]+a3[0:30])
        features(up[-30:]+down[0:30])
        features(a3[-30:]+down[0:30])
        #6
        features(up[-50:]+a3[0:50])
        features(up[-50:]+down[0:50])
        features(a3[-50:]+down[0:50])
        #7
        features(up[-20:])
        features(up[-10:])
        features(down[0:10])
        features(down[0:20])
        print("0")




 #!/usr/bin/python
# -*- coding: utf-8 -*-
#coding=utf-8

#from Bio import SeqIO
#
#seqs={}
#gene={}
##for seq_record in SeqIO.parse("Araport11_cdna_onlyname.fasta", "fasta"):
#for seq_record in SeqIO.parse("example.fasta", "fasta"):
#	seqs[seq_record.id]=seq_record.seq._data
import re
import sys

sys.setrecursionlimit(900)
#定义k-mer长度
#k=25
#定义一个字典存图
#定义一个字典存bubble



#取最小重复序列，即k
def mink(seq,k,flag,pos):
	pos = pos
	if flag == 1:
		for i in range(pos,len(seq)-k+1):
			if seq.count(seq[i:i+k]) > 1 or seq[i:i+k] in seq[i+1:i+k+k]:
				flag=1
				pos = i
				break
			else:
				flag=0
		k = mink(seq,k+1,flag,pos)
		return k
	else:
		return k

#构建图，把两个序列按顺序，放到nodes字典中
def construteDG(seq1, seq2):
	#if len(seq1) < len(seq2):
	#	seq1,seq2 = seq2,seq1
	nodes={}
	bubbles={}
	k=10
	flag = 1
	k=max(mink(seq1.upper(),k,flag,0),mink(seq2.upper(),k,flag,0))
	#k=mink(seq1,k,flag)+1
	tail="Q"*(k-1)+"PPPP"
	head="Q"*(k)
	seq1 = head+seq1.upper()+tail
	seq2 = head+seq2.upper()+tail
	#print(seq1,seq2,sep="\n")
	#该循环构建第一个序列
	for i in range(0,int(f'{len(seq1)}')-k+1):
		#将第一个节点的父节点，设为空值
		if i==0:
			parent = "null"
		else:
			parent = seq1[i-1:i-1+k]
		
		#将最后一个节点的子节点，设为空值
		if i== int(f'{len(seq1)}')-k:
			son = "null"
		else:
			son = seq1[i+1:i+1+k]
		
		node=seq1[i:i+k]

		if node in nodes:
			print("存在重复kmer"+node+"，请调大k值"+str(k)+","+str(len(nodes)))
			print(seq1)
			print(seq2)
			break
		else:
			nodes[node]={
				'source':'a',
				'seq':node,
				'length':k,
				'position':i+1,#以node的末位置来表示在序列中的位置，帽子的位置都认为是负数
				'parent':parent,
				'son': son,}
		#print(nodes[node])
	#在不在之前的node里面都是一个分岔点，是一个新路的起止点
	path_st=path_cv=0
	flag="st"
	bu_number=1	
	#该循环构建第二个序列
	for i in range(0,int(f'{len(seq2)}')-k):
		#将第一个节点的父节点，设为空值
		if i==0:
			parent = "null"
		else:
			parent = seq2[i-1:i-1+k]
		
		#将最后一个节点的子节点，设为空值
		if i== int(f'{len(seq2)}')-k:
			son = "null"
		else:
			son = seq2[i+1:i+1+k]
		node=seq2[i:i+k]
		#print(node)
		#判断来源,是否在序列1中已经存在
		if node in nodes:
			#print(nodes[node])
			if flag=="cv":#一个bubble的结束
				flag="st"
				#一个bubble的所有信息
				#length1是这个node重新出现的地方，也就是一个bubble的结束位置，减去，记录在bubbles字典里的这个bubble的的起始点，再减去，k
				#
				length1 = nodes[node]["position"] - bubbles["bubble"+str(bu_number)]["start_pos1"]-k
				#
				length2 = i - bubbles["bubble"+str(bu_number)]["start_pos2"] - k + 1
				#
				end_pos1= nodes[node]["position"]-k
				#
				end_pos2= i+1-k
				phase=0
				#序列中总会存在bubble的前头和尾巴有相同序列的 问题，所以在bubble结束时，总会出现位点前移的问题，但是这部分序列已经算在了开头，所以要进行矫正
				#如果是第一个序列的重复，那么length1会小于0，这时需要对两个end_pos都减去这个负值，即加上这个重复片段的长度，此时bubble的sink点会移到正确位置
				if length1 < 0 :
					#该变量是为了输出查看
					phase=length1
					#将该bubble的sink位点，后移至正确位置
					end_pos1 = end_pos1 - length1
					end_pos2 = end_pos2 - length1
					length2 = length2 - length1
					length1 = 0
					#矫正完之后该路径为0
				if length2 < 0 :
					phase=length2
					end_pos1 = end_pos1 - length2
					end_pos2 = end_pos2 - length2
					length1 = length1 - length2
					length2 = 0
				if abs(phase) > k:
					continue
				bubbles["bubble"+str(bu_number)]["end_pos1"]=end_pos1
				bubbles["bubble"+str(bu_number)]["end_pos2"]=end_pos2
				bubbles["bubble"+str(bu_number)]["phase"]=phase
				bubbles["bubble"+str(bu_number)]["end_node"]=node
				bubbles["bubble"+str(bu_number)]["length1"]= length1
				bubbles["bubble"+str(bu_number)]["length2"]= length2
				
				path_st=0
				bu_number+=1
			else:#一个桥的延续
				path_st+=1
				flag="st"
		else:
			if flag=="st":#一个桥的结束，一个bubble的开始
				flag="cv"
				#bubble的前一个分岔点的信息
				bubbles["bubble"+str(bu_number)]={
				"start_node":seq2[i-1:i-1+k],
				"start_pos1":nodes[seq2[i-1:i-1+k]]["position"],
				"start_pos2":i
				}
				path_cv=0
			else:#一个bubble的延续
				path_cv+=1
				flag="cv"
	#print(bubbles)
	del nodes
	return bubbles
'''
截取特征片段
	for x in bubbles:
		if bubbles[x]["length2"] == 0:
			print("AA",end=",")
			print(seq1[bubbles[x]["start_pos1"]-30:bubbles[x]["start_pos1"]],end=",")
			print(seq1[bubbles[x]["start_pos1"]:bubbles[x]["end_pos1"]],end=",")
			print(seq1[bubbles[x]["end_pos1"]:bubbles[x]["end_pos1"]+30],end="\n")
		elif bubbles[x]["length1"] == 0:
			print("AA",end=",")
			print(seq2[bubbles[x]["start_pos2"]-30:bubbles[x]["start_pos2"]],end=",")
			print(seq2[bubbles[x]["start_pos2"]:bubbles[x]["end_pos2"]],end=",")
			print(seq2[bubbles[x]["end_pos2"]:bubbles[x]["end_pos2"]+30],end="\n")
'''
	

#判断bubble的性质
def justify(bubbles,seq1,seq2):
	answer=""
	if seq1+seq2 in ass:
		#print(seq1,seq2,sep="\t")
		answer = ass[seq1+seq2]
	#print(answer)
	num=0
	snpnum=0
	#print(bubbles)
	for bubble in bubbles:
		bubble=bubbles[bubble]
		#print(bubble)
		if "length1" in bubble:
			if bubble["length1"] and bubble["length1"] == 1 and bubble["length2"] == 1:
				#print("this is a snp")
				snpnum=snpnum+1
			else:
				if str(bubble['start_pos1']) in answer or str(bubble['start_pos2']) in answer or str(bubble['end_pos1']) in answer or str(bubble['end_pos2']) in answer:
					#print(bubble,answer,sep="\t")
					num=num+1
	#print(num)
	return num,snpnum

#ass = {}
def getnewsuppa():
	ass = {}
#	suppa = open("./newall.ioe","r")
	suppa = open("./changesupparesult.out","r")
	
	while 1:
		lines = suppa.readlines(100000)
		if not lines:
			break
		for line in lines:
			line=line.strip('\n')
			lines=line.split('\t')
			bothname=""
			if lines[2] < lines[3]:
				bothname = lines[2]+lines[3]
			else:
				bothname = lines[3]+lines[2]
			if bothname in ass:
				ass[bothname].append([int(lines[4]),int(lines[5]),int(lines[6]),int(lines[7]),lines[0]])
			else:
				ass[bothname]=[[int(lines[4]),int(lines[5]),int(lines[6]),int(lines[7]),lines[0]]]
	return ass

def getnewgff():
	gffs = {}
	gff = open("./new.gff","r")
	tmp = "a"

	while 1:
		lines = gff.readlines(100000)
		if not lines:
			break
		for line in lines:
			hangs=line.split('\t')
			#searchObj = re.search( r'^AT', line, re.M|re.I)
			if hangs[0][0:2] == "AT":
				tmp = hangs[0]
				gffs[tmp] = [hangs[1]]
			else:
				gffs[tmp].append(hangs[0])
				gffs[tmp].append(hangs[1])
	return gffs

def gen2trans():#基因组和转录组的位置互换	transpos = {}
	gff = open("./new.gff","r")
	tmp = "a"
	transpos["a"]="开始"

	while 1:
		lines = gff.readlines(100000)
		if not lines:
			break
		for line in lines:
			line=line.strip('\n')
			hangs=line.split('\t')
			#searchObj = re.search( r'^AT', line, re.M|re.I)
			if hangs[0][0:2] == "AT":
				#print(transpos[tmp])
				tmp = hangs[0]#转录本名
				transpos[tmp] = [hangs[1]]
			else:
				transpos[tmp+hangs[2]] = hangs[0]#转录本名+基因组上的位置，等于转录组上的位置
				transpos[tmp+hangs[3]] = hangs[1]
	return transpos


def changposforvista(name,pos):
	tr_poss=[]
	#print(name,pos)
	for x in pos.split(","):
		if name+str(x) in transpos:
			tr_poss.append(transpos[name+str(x)])
			#print(transpos[name+str(x)])
	#print(tr_poss)
	return(tr_poss)
def vista():
    vista_num=0
    bubble_num=0
    ioe = open("./astalavista.gtf","r")
    vistass={}
    while 1:
        lines = ioe.readlines(100000)
        if not lines:
            break
        for line in lines:
            vista_num =vista_num +1
            line=line.strip('\n')
            #as_code "1^2- , 0"; transcript1_id "AT1G36180.2"; splice_chain1 ""; transcript2_id "AT1G36180.3"; splice_chain2 "13547775,13547864"; 
            searchObj = re.search( r'.*?transcript1_id "(.*?)"; splice_chain1 "(.*?)"; transcript2_id "(.*?)"; splice_chain2 "(.*?)";.*', line, re.M|re.I)
            if searchObj:
                #print(searchObj.group(0))
                tran1=searchObj.group(1)
                tran2=searchObj.group(3)
                bothname=""
                if tran1 < tran2:
                    bothname = tran1 + tran2
                else:
                    bothname = tran2 + tran1
                if bothname in vistass:
                	if searchObj.group(2) == "":
                		vistass[bothname].extend(changposforvista(searchObj.group(3),searchObj.group(4)))
                	elif searchObj.group(4) == "":
                		vistass[bothname].extend(changposforvista(searchObj.group(1),searchObj.group(2)))
                	else:
                		vistass[bothname].extend(changposforvista(searchObj.group(1),searchObj.group(2)))
                		vistass[bothname].extend(changposforvista(searchObj.group(3),searchObj.group(4)))
                else:
                	if searchObj.group(2) == "":
                		vistass[bothname]=changposforvista(searchObj.group(3),searchObj.group(4))
                	elif searchObj.group(4) == "":
                		vistass[bothname]=changposforvista(searchObj.group(1),searchObj.group(2))
                	else:
                		vistass[bothname]=changposforvista(searchObj.group(1),searchObj.group(2))
                		vistass[bothname].extend(changposforvista(searchObj.group(3),searchObj.group(4)))
    return vistass
#transpos=gen2trans()
#print(transpos)

				
def mink2(seq,k,flag,pos):
	#print(seq)
	pos = pos
	if flag == 1:
		for i in range(pos,len(seq)-k+1):
			if seq.count(seq[i:i+k]) > 1 or seq[i:i+k] in seq[i+1:i+k+k]:
				flag=1
				pos = i
				break
			else:
				flag=0
			#print(i,k,seq[i:i+k],seq.count(seq[i:i+k]))
		k = mink2(seq,k+1,flag,pos)
		return k
	else:
		return k
       
			


def secondconstruteDG(start1,start2,bubble_number,seq1, seq2):
	#if len(seq1) < len(seq2):
	#	seq1,seq2=seq2,seq1
	#print(bubble_number)
	#print(seq1)
	#print('')
	#print(seq2)
	nodes={}
	bubbles={}

	#print(type(seq1),type(seq2))
	k=4
	flag = 1
	#print("ca",seq1, seq2)
	k=max(mink2(seq1,k,flag,0),mink2(seq2,k,flag,0))
	#k=mink2(seq1,k,flag,0)-1
	#print(k)
	tail="Q"*(k-1)+"PPPP"
	head="Q"*(k)
	seq1 = head+seq1+tail
	seq2 = head+seq2+tail
	#print(seq1,seq2,sep="\n")
	#该循环构建第一个序列
	for i in range(0,int(f'{len(seq1)}')-k+1):
		#将第一个节点的父节点，设为空值
		if i==0:
			parent = "null"
		else:
			parent = seq1[i-1:i-1+k]
		
		#将最后一个节点的子节点，设为空值
		if i== int(f'{len(seq1)}')-k:
			son = "null"
		else:
			son = seq1[i+1:i+1+k]
		
		node=seq1[i:i+k]

		if node in nodes:
			print("存在重复kmer"+node+"，请调大k值"+str(k)+str(len(nodes)))
			print(seq1)
			break
		else:
			nodes[node]={
				'source':'a',
				'seq':node,
				'length':k,
				'position':i+1,#以node的末位置来表示在序列中的位置，帽子的位置都认为是负数
				'parent':parent,
				'son': son,}
		#print(nodes[node])
	#在不在之前的node里面都是一个分岔点，是一个新路的起止点
	path_st=path_cv=0
	flag="st"
	bu_number=1
	#该循环构建第二个序列
	for i in range(0,int(f'{len(seq2)}')-k):
		#将第一个节点的父节点，设为空值
		if i==0:
			parent = "null"
		else:
			parent = seq2[i-1:i-1+k]
		
		#将最后一个节点的子节点，设为空值
		if i== int(f'{len(seq2)}')-k:
			son = "null"
		else:
			son = seq2[i+1:i+1+k]
		node=seq2[i:i+k]
		#print(node)
		#判断来源,是否在序列1中已经存在
		if node in nodes:
			#print(nodes[node])
			if flag=="cv":#一个bubble的结束
				flag="st"
				#一个bubble的所有信息
				#length1是这个node重新出现的地方，也就是一个bubble的结束位置，减去，记录在bubbles字典里的这个bubble的的起始点，再减去，k
				#
				length1 = nodes[node]["position"] - bubbles[bubble_number+"."+str(bu_number)]["start_pos1"]-k
				#
				length2 = i - bubbles[bubble_number+"."+str(bu_number)]["start_pos2"] - k + 1
				#
				end_pos1= nodes[node]["position"]-k
				#
				end_pos2= i+1-k
				phase=0
				#序列中总会存在bubble的前头和尾巴有相同序列的问题，所以在bubble结束时，总会出现位点前移的问题，但是这部分序列已经算在了开头，所以要进行矫正
				#如果是第一个序列的重复，那么length1会小于0，这时需要对两个end_pos都减去这个负值，即加上这个重复片段的长度，此时bubble的sink点会移到正确位置
				if length1 < 0 :
					#该变量是为了输出查看
					phase=length1
					#将该bubble的sink位点，后移至正确位置
					end_pos1 = end_pos1 - length1
					end_pos2 = end_pos2 - length1
					length2 = length2 - length1
					length1 = 0
					#矫正完之后该路径为0
				if length2 < 0 :
					phase=length2
					end_pos1 = end_pos1 - length2
					end_pos2 = end_pos2 - length2
					length1 = length1 - length2
					length2 = 0
				bubbles[bubble_number+"."+str(bu_number)]["end_pos1"]=end_pos1
				bubbles[bubble_number+"."+str(bu_number)]["end_pos2"]=end_pos2
				bubbles[bubble_number+"."+str(bu_number)]["phase"]=phase
				bubbles[bubble_number+"."+str(bu_number)]["end_node"]=node
				bubbles[bubble_number+"."+str(bu_number)]["length1"]= length1
				bubbles[bubble_number+"."+str(bu_number)]["length2"]= length2
				
				path_st=0
				bu_number+=1
			else:#一个桥的延续
				path_st+=1
				flag="st"
		else:
			if flag=="st":#一个桥的结束，一个bubble的开始
				flag="cv"
				#bubble的前一个分岔点的信息
				bubbles[bubble_number+"."+str(bu_number)]={
				"start_node":seq2[i-1:i-1+k],
				"start_pos1":nodes[seq2[i-1:i-1+k]]["position"],
				"start_pos2":i
				}
				path_cv=0
			else:#一个bubble的延续
				path_cv+=1
				flag="cv"
	#print(start1)
	#print(start2)
	for x in bubbles:
		bubbles[x]["start_pos1"]=int(bubbles[x]["start_pos1"])+int(start1)-1
		bubbles[x]["start_pos2"]=int(bubbles[x]["start_pos2"])+int(start2)-1
		bubbles[x]["end_pos1"]=int(bubbles[x]["end_pos1"])+int(start1)-1
		bubbles[x]["end_pos2"]=int(bubbles[x]["end_pos2"])+int(start2)-1
	#print(bubbles)
	del nodes
	return bubbles

'''

def getnewsuppa():
	ass = {}
#	suppa = open("./newall.ioe","r")
	suppa = open("./changesupparesult.out","r")
	
	while 1:
		lines = suppa.readlines(100000)
		if not lines:
			break
		for line in lines:
			lines=line.split('\t')
			if lines[2]<lines[3]:
				if lines[2]+lines[3] in ass:
					ass[lines[2]+lines[3]].extend([lines[4],lines[5],lines[6],lines[7],lines[0]])
				else:
					ass[lines[2]+lines[3]]=[lines[4],lines[5],lines[6],lines[7],lines[0]]
			else:
				if lines[3]+lines[2] in ass:
					ass[lines[3]+lines[2]].extend([lines[4],lines[5],lines[6],lines[7],lines[0]])
				else:
					ass[lines[3]+lines[2]]=[lines[4],lines[5],lines[6],lines[7],lines[0]]
	return ass

'''

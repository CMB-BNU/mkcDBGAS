#!/user/local/bin/python
#coding=utf-8

#引入包
import re
import sys
import numpy as np
from sys import argv
import pandas as pd
import joblib
import datetime

now_time = datetime.datetime.now()
print(now_time)
print(argv[1])

#加载训练数据
data_path= argv[1] #数据文件的路径
data = np.loadtxt(data_path,                                #数据文件路径
                  dtype=float,                              #数据类型
                  delimiter=',')                            #数据分隔符

#print(data.shape)
#head_path='feature.head'  #数据文件的路径
#head = np.loadtxt(head_path,dtype=str,delimiter=',')
#astype,feature = np.split(head,(1,),axis=1)
#数据分割
y, x = np.split(data,                                       #要切分的数组
                (0,),                                       #沿轴切分的位置，第5列开始往后为y
                axis=1)                                     #代表纵向分割，按列分割

#取出模型
scaler=joblib.load("./models/"+argv[2]+"_scaler.pkl")
featureselector=joblib.load("./models/"+argv[2]+"_featureselector.pkl")
clf=joblib.load("./models/"+argv[2]+"_clf.pkl")
#标准化
x = scaler.transform(x)
#特征选择
x = featureselector.transform(x)


y_pre = clf.predict(x)
y_pre = y_pre[:,np.newaxis]
y_pro = clf.predict_proba(x)
y_ppp=[]
for s in y_pre:
    if s[0] == 0:
        y_ppp.append(["AA"])
    if s[0] == 1:
        y_ppp.append(["AD"])
    if s[0] == 2:
        y_ppp.append(["ES"])
    if s[0] == 3:
        y_ppp.append(["IR"])
y_pre = np.array(y_ppp)

event_path=argv[3]  #数据文件的路径
event = np.loadtxt(event_path,dtype=str,delimiter=',')


c = np.hstack((event,y_pre,y_pro))
d = pd.DataFrame(c)
d.to_csv(argv[4],index=False,sep=',')

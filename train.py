#!/user/local/bin/python
#coding=utf-8

#引入包
import re
import sys
import numpy as np
from sys import argv
import pandas as pd
from matplotlib import colors
from sklearn import svm
from sklearn.svm import SVC
from sklearn import model_selection
import matplotlib as mpl
from sklearn.metrics import confusion_matrix, precision_score, accuracy_score,recall_score, f1_score,roc_auc_score
import joblib
import pickle
from sklearn.metrics import classification_report
import datetime
from sklearn import preprocessing
from xgboost import XGBClassifier
from sklearn.ensemble import ExtraTreesClassifier
from imblearn.over_sampling import BorderlineSMOTE
from sklearn.feature_selection import SelectFromModel
#*************将字符串转为整型，便于数据加载***********************
def as_type(s):
    it = {b'A3':0, b'A5':1, b'SE':2, b'RI':3}
    return it[s]

now_time = datetime.datetime.now()
print(now_time)
print(argv[1])

#加载训练数据
data_path='/panfs4/home/zhangqb/mkcDBG/feature/seqara50.feature1-7'  #数据文件的路径
#data_path='/panfs4/home/zhangqb/modeltrain/human/seqhuman50.feature1-7'
data = np.loadtxt(data_path,                                #数据文件路径
                  dtype=float,                              #数据类型
                  delimiter=',',                            #数据分隔符
                  converters={0:as_type})                 #将第5列使用函数as_type进行转换

#print(data.shape)
#数据分割
y, x = np.split(data,                                       #要切分的数组
                (1,),                                       #沿轴切分的位置，第5列开始往后为y
                axis=1)                                     #代表纵向分割，按列分割
#x = x[:, 0:2]                                               #在X中我们取前两列作为特征，为了后面的可视化。x[:,0:4]代表第一维(行)全取，第二维(列)取0~2
print(x.shape)
print(y)

#标准化
scaler = preprocessing.StandardScaler().fit(x)
x = scaler.transform(x)

#保存
joblib.dump(scaler, "/panfs4/home/zhangqb/mkcDBG/models/"+argv[1]+"_scaler.pkl")

#特征选择
featureselector = SelectFromModel(estimator=ExtraTreesClassifier()).fit(x, y.ravel())
x = featureselector.transform(x)
featureimportance = featureselector.estimator_.feature_importances_
print(featureimportance.shape)
featureimportance = featureimportance.reshape(1,-1)
#print(featureimportance)
#过采样
cc = BorderlineSMOTE(random_state=0,kind="borderline-2")
x, y = cc.fit_resample(x, y)


featureimportance = featureselector.transform(featureimportance)
#保存
joblib.dump(featureselector, "/panfs4/home/zhangqb/mkcDBG/models/"+argv[1]+"_featureselector.pkl")
print(x.shape)
x_train,x_test,y_train,y_test=model_selection.train_test_split(x,              #所要划分的样本特征集
                                                               y,              #所要划分的样本结果
                                                               random_state=1, #随机数种子
                                                               test_size=0.3)  #测试样本占比
#print(x_train.shape)
# 2.定义模型：SVM模型定义
clf = XGBClassifier()
#***********************训练模型*****************************
def train(clf,x_train,y_train):
    clf.fit(x_train,         #训练集特征向量
            y_train.ravel()) #训练集目标值

# 3.训练SVM模型
train(clf,x_train,y_train)
#保存
joblib.dump(clf, "/panfs4/home/zhangqb/mkcDBG/models/"+argv[1]+"_clf.pkl")
print("train done")


def show_accuracy(a, b, tip):
    acc = a.ravel() == b.ravel()
    print('%s Accuracy:%.3f' %(tip, np.mean(acc)))

def print_accuracy(clf,x_train,y_train,x_test,y_test):
    #分别打印训练集和测试集的准确率  score(x_train,y_train):表示输出x_train,y_train在模型上的准确率
    print('trianing prediction:%.3f' %(clf.score(x_train, y_train)))
    print('test data prediction:%.3f' %(clf.score(x_test, y_test)))
    #原始结果与预测结果进行对比   predict()表示对x_train样本进行预测，返回样本类别
    #show_accuracy(clf.predict(x_train), y_train, 'traing data')
    #show_accuracy(clf.predict(x_test), y_test, 'testing data')
    print(classification_report(y_test,clf.predict(x_test),digits=5))
    #计算决策函数的值，表示x到各分割平面的距离
    #print('decision_function:\n', clf.decision_function(x_train))

print_accuracy(clf,x_train,y_train,x_test,y_test)

# 4.获取特征
head_path='/panfs4/home/zhangqb/mkcDBG/feature/seqara50.feature1-7.head'  #数据文件的路径
head = np.loadtxt(head_path,dtype=str,delimiter=',')
astype,feature = np.split(head,(1,),axis=1)
print(feature.shape)
feature = featureselector.transform(feature)
print(feature[0])
print(featureimportance)
now_time = datetime.datetime.now()
print(now_time)

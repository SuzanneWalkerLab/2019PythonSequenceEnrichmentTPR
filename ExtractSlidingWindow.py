#!/usr/bin/env python2
#sliding window enrichment analysis of HeLa extract data

import numpy as np
import pandas as pd
from Functions.py import AAwindow
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt

#read in data
data = pd.read_csv('data.csv', header = 0)

#find the score for each peptide
Cols_of_interest = []
for i in ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','W','Y']:
    for j in range(16):
        Col = i+'_C'+str(j)+'-'+str(j+2)
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = (j,j+3), aa = i),
                           axis=1)

#generate AHNAK-free data
dataNA = data.loc[data.Accession != 'Q09666',].reset_index(drop = True)

#analyze enrichment
WTonly = (dataNA.WT > 0) & (dataNA.Mutant == 0)
D2Aonly =  (dataNA.Mutant > 0) & (dataNA.WT == 0)
dataList = [dataNA,#0
       dataNA.loc[WTonly,:],#1 WTonly
       dataNA.loc[D2Aonly,:]#2 D2Aonly
       ]
comparisons = [(1,2)]
compDesc = ['WTonly_v_D2Aonly']
labels = ['All','WTonly','D2Aonly']
mannU = mannUpval(Cols_of_interest, dataList, comparisons, compDesc, labels, mean = True, median = False,
              fold = True, diff = True, storeU = False)

#save "hits" from 3aa analysis
where = (mannU.WTonly_v_D2Aonly_pval < 0.05) &
          (mannU.WTonly_v_D2Aonly_fold.abs() > 1.5)
mannU.loc[where,].to_csv('3aaHits.csv',index = False)

#plot volcano plot

ax = plt.gca()
mannU.loc[where & mannU.Cols.str.count('K_C'),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'blue',
                         ax = ax)
mannU.loc[where &mannU.Cols.str.count('R_C'),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'blue',
                         ax = ax)
mannU.loc[where&(mannU.Cols.str.count('V_C') | mannU.Cols.str.count('L_C')),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'black',
                         ax = ax)
mannU.loc[where&mannU.Cols.str.findall('S_C'),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'green',
                         ax = ax)
mannU.loc[where&mannU.Cols.str.findall('T_C'),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'green',
                         ax = ax)
mannU.loc[where&mannU.Cols.str.findall('A_C'),:].plot(
    x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',kind = 'scatter',color = 'pink',
                         ax = ax)
mannU.loc[~where,:].plot(x='WTonly_v_D2Aonly_log2fold',y='pWTonly_v_D2Aonly',
                           color = 'gray',kind = 'scatter', ax = ax)



plt.title('Comparison of WT specific versus D2A specific sites')
plt.xlabel('Log2 Fold Change')
plt.yticks([0,1,2,3],[1,0.1,0.01,0.001])
plt.ylabel('p-value')
#plt.legend(['K','R','V','S','T','A','nonsignificant']) #covers significant proportion, include to follow coloring.
#colors adjusted post-hoc in illustrator to match color scheme
plt.savefig('./WTonly_v_D2Aonly_3s.pdf')

#expanded window analysis:
#increase or decrease window size, combine or don't
Cols_of_interest = []
#Ser/Thr early
for i in ['S','T','ST']:
    for j in [(5,6),(5,7),(5,8),(5,9),
              (6,7),(6,8),(6,9),
             (7,8),(7,9),
             (8,9)]:
        Col = i+'_C'+str(j[0])+'-'+str(j[1]-1)
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = j, aa = i),
                           axis=1)
        
#Ser/Thr late
for i in ['S','T','ST']:
    for j in [(11,12),(11,13),(11,14),(11,15),
              (12,13),(12,14),(12,15),
             (13,14),(13,15),
             (14,15)]:
        Col = i+'_C'+str(j[0])+'-'+str(j[1])
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = j, aa = i),
                           axis=1)
        
#Ala
for i in ['A']:
    for j in [(11,12),(11,13),(11,14),
              (12,13),(12,14),
             (13,14)]:
        Col = i+'_C'+str(j[0])+'-'+str(j[1])
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = j, aa = i),
                           axis=1)

#Arg/Lys
for i in ['R','K','RK']:
    for j in [(7,8),(7,9),(7,10),(7,11),(7,12),(7,13),(7,14),
              (8,9),(8,10),(8,11),(8,12),(8,13),(8,14),
             (9,10),(9,11),(9,12),(9,13),(9,14),
             (10,11),(10,12),(10,13),(10,14),
             (11,12),(11,13),(11,14),
             (12,13),(12,14),
             (13,14)]:
        Col = i+'_C'+str(j[0])+'-'+str(j[1]-1)
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = j, aa = i),
                           axis=1)

#Val/Ile/Leu 2-5, 8-11
for i in ['V','VI','VL','VIL']:
    for j in [(2,3),(2,4),(2,5),
              (3,4),(3,5),
             (4,5),
             (8,9),(8,10),(8,11),
             (9,10),(9,11),
             (10,11)]:
        Col = i+'_C'+str(j[0])+'-'+str(j[1])
        Cols_of_interest += [Col]
        data[Col] = data.apply(
            lambda row: AAwindow(row['Sequence'], row['Site'],Window = j, aa = i),
                           axis=1)

#re-generate AHNAK-free data
dataNA = data.loc[data.Accession != 'Q09666',].reset_index(drop = True)



data = [dataNA.loc[WTonly,:],#0
        dataNA.loc[D2Aonly,:]]#1] 
comparisons = [(0,1)]#13
compDesc = ['WTonly_v_D2Aonly']
labels = ['WTonly','D2Aonly']
mannU_Subsetted = mannUpval(Cols_of_interest, data, comparisons, compDesc, labels, mean = True, median = False,
              fold = True, diff = True, storeU = False)
#table output of expanded windows
mannU_Subsetted.loc[mannU_Subsetted.WTonly_v_D2Aonly_pval < 0.001,].to_csv('FinalWindows.csv',index = False)
#bar charts of top windows

#what are values?
col = 'RK_C7-11'
colname = 'RK 7-11'
WTonly = (data.Mutant == 0) & (data.WT > 0)
denom = sum(WTonly)
sites = WTonly
name = 'Glycosylated by WTonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'
print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname
print float(sum(data.loc[sites,:][col]==3))/denom*100,'% 3 matches '+colname
print float(sum(data.loc[sites,:][col]==4))/denom*100,'% 4 matches '+colname
print float(sum(data.loc[sites,:][col]==5))/denom*100,'% 5 matches '+colname

D2Aonly = (data.Mutant > 0) & (data.WT == 0)
denom = sum(D2Aonly)
sites = D2Aonly
name = 'Glycosylated by D2Aonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'

print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname
print float(sum(data.loc[sites,:][col]==3))/denom*100,'% 3 matches '+colname
print float(sum(data.loc[sites,:][col]==4))/denom*100,'% 4 matches '+colname
print float(sum(data.loc[sites,:][col]==5))/denom*100,'% 5 matches '+colname


col = 'ST_C6-8'
colname = 'ST 6-8'
denom = sum(WTonly)
sites = WTonly
name = 'Glycosylated by WTonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'
print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname
print float(sum(data.loc[sites,:][col]==3))/denom*100,'% 3 matches '+colname


denom = sum(D2Aonly)
sites = D2Aonly
name = 'Glycosylated by D2Aonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'

print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname
print float(sum(data.loc[sites,:][col]==3))/denom*100,'% 3 matches '+colname


col = 'S_C13-14'
colname = 'S 13-14'
denom = sum(WTonly)
sites = WTonly
name = 'Glycosylated by WTonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'
print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname


denom = sum(D2Aonly)
sites = D2Aonly
name = 'Glycosylated by D2Aonly'
print "______All " + name + " sites from glycosite mapping______"
print denom, ' ' + name + ' Sites'

print float(sum(data.loc[sites,:][col]==1))/denom*100,'% 1 matches '+colname
print float(sum(data.loc[sites,:][col]==2))/denom*100,'% 2 matches '+colname

#make bar charts based upon above numbers (note re-colored in illustrator post-creation)


#bar plots stacked

RK_1 = np.array([45,28.1])
RK_2 = np.array( [9.9,6.7])
RK_3 = np.array( [4.5,0])





N = 2
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars

p1 = plt.bar(ind, RK_1, width)
p2 = plt.bar(ind, RK_2, width,bottom=RK_1)
p3 = plt.bar(ind, RK_3, width,bottom=RK_1+RK_2)




plt.ylabel('Percent Glycosites')
plt.title('Distribution of Pattern Matching for Extract-Glycosylated Proteins (RK 7-11)')
plt.xticks(ind, ('WT only',  'D2A only'))
plt.ylim((0,70))

plt.legend((p1[0], p2[0], p3[0]),
           ('1 Arg/Lys', '2 Arg/Lys', '3 Arg/Lys'))

plt.savefig('RK7-11BarGraph_withAHNAK.pdf')
plt.clf()

#ST6-8
ST_1 = np.array([23.4,43.7])
ST_2 = np.array( [5.4,12.6])
ST_3 = np.array( [0,0.74])


N = 2
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars

p1 = plt.bar(ind, ST_1, width)
p2 = plt.bar(ind, ST_2, width,bottom=ST_1)
p3 = plt.bar(ind, ST_3, width,bottom=ST_1+ST_2)




plt.ylabel('Percent Glycosites')
plt.title('Distribution of Pattern Matching for Extract-Glycosylated Proteins (ST 6-8)')
plt.xticks(ind, ('WT only',  'D2A only'))
plt.ylim((0,60))

plt.legend((p1[0], p2[0], p3[0]),
           ('1 Ser/Thr', '2 Ser/Thr', '3 Ser/Thr'))

plt.savefig('ST6-8BarGraph_withAHNAK.pdf')
plt.clf()

S_1 = np.array([6.7,23.7])
S_2 = np.array( [0,2.3])





N = 2
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars

p1 = plt.bar(ind, S_1, width)
p2 = plt.bar(ind, S_2, width,bottom=S_1)




plt.ylabel('Percent Glycosites')
plt.title('Distribution of Pattern Matching for Ser/Thr from Extract-Glycosylated Proteins (S 13-14)')
plt.xticks(ind, ('WT only',   'D2A only'))
plt.ylim((0,30))

plt.legend((p1[0], p2[0]),
           ('1 Ser', '2 Ser'))

plt.savefig('S13-14BarGraph.pdf')

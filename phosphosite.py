#!/usr/bin/env python2
#sliding window enrichment analysis of phosphosite data

import numpy as np
import pandas as pd
import requests
import re
from Functions.py import GetSeqUniprot
from Functions.py import AAwindow
from Functions.py import findST
from scipy.stats import mannwhitneyu
from matplotlib import pyplot as plt
from statsmodels.stats.multitest import multipletests
from Functions.py import plotLines


#read in data; downloaded from PhosphoSitePlus on Aug 8, 2018
phosphosite = pd.read_table('O-GlcNAc_site_dataset', header = 2)
BASE = "http://www.uniprot.org"
KB_ENDPOINT = "/uniprot/"


phosphosite['Prot_Seq'] = phosphosite.ACC_ID.apply(GetSeqUniprot)


phosphosite['Site'] = phosphosite.MOD_RSD.apply(lambda x: int(re.search('[ST]([0-9]+)\-gl',x).group(1)))
#to replace seq for obsolete sequence
for x in range(1201,1213):
    phosphosite.loc[x,'Prot_Seq']= 'MADPAAGPPPSEGEESTVRFARKGALRQKNVHEVKNHKFTARFFKQPTFCSHCTDFIWGFGKQGFQCQVCCFVVHKRCHEFVTFSCPGADKGPASDDPRSKHKFKIHTYSSPTFCDHCGSLLYGLIHQGMKCDTCMMNVHKRCVMNVPSLCGTDHTERRGRIYIQAHIDREVLIVVVRDAKNLVPMDPNGLSDPYVKLKLIPDPKSESKQKTKTIKCSLNPEWNETFRFQLKESDKDRRLSVEIWDWDLTSRNDFMGSLSFGISELQKAGVDGWFKLLSQEEGEYFNVPVPPEGSEGNEELRQKFERAKIGQGTKAPEEKTANTISKFDNNGNRDRMKLTDFNFLMVLGKGSFGKVMLSERKGTDELYAVKILKKDVVIQDDDVECTMVEKRVLALPGKPPFLTQLHSCFQTMDRLYFVMEYVNGGDLMYHIQQVGRFKEPHAVFYAAEIAIGLFFLQSKGIIYRDLKLDNVMLDSEGHIKIADFGMCKENIWDGVTTKTFCGTPDYIAPEIIAYQPYGKSVDWWAFGVLLYEMLAGQAPFEGEDEDELFQSIMEHNVAYPKSMSKEAVAICKGLMTKHPGKRLGCGPEGERDIKEHAFFRYIDWEKLERKEIQPPYKPKARDKRDTSNFDKEFTRQPVELTPTDKLFIMNLDQNEFAGFSYTNPEFVINVMADPAAGPPPSEGEESTVRFARKGALRQKNVHEVKNHKFTARFFKQPTFCSHCTDFIWGFGKQGFQCQVCCFVVHKRCHEFVTFSCPGADKGPASDDPRSKHKFKIHTYSSPTFCDHCGSLLYGLIHQGMKCDTCMMNVHKRCVMNVPSLCGTDHTERRGRIYIQAHIDREVLIVVVRDAKNLVPMDPNGLSDPYVKLKLIPDPKSESKQKTKTIKCSLNPEWNETFRFQLKESDKDRRLSVEIWDWDLTSRNDFMGSLSFGISELQKAGVDGWFKLLSQEEGEYFNVPVPPEGSEGNEELRQKFERAKIGQGTKAPEEKTANTISKFDNNGNRDRMKLTDFNFLMVLGKGSFGKVMLSERKGTDELYAVKILKKDVVIQDDDVECTMVEKRVLALPGKPPFLTQLHSCFQTMDRLYFVMEYVNGGDLMYHIQQVGRFKEPHAVFYAAEIAIGLFFLQSKGIIYRDLKLDNVMLDSEGHIKIADFGMCKENIWDGVTTKTFCGTPDYIAPEIIAYQPYGKSVDWWAFGVLLYEMLAGQAPFEGEDEDELFQSIMEHNVAYPKSMSKEAVAICKGLMTKHPGKRLGCGPEGERDIKEHAFFRYIDWEKLERKEIQPPYKPKARDKRDTSNFDKEFTRQPVELTPTDKLFIMNLDQNEFAGFSYTNPEFVINV'


#generate background list of all Ser/Thr



Sequences = []
Accessions = []
FakeSite = []
RK = []
DE = []
for seq in data.Prot_Seq.unique():
    for res in findST(seq):
        Sequences += [seq]
        Accessions += [getAccession(data,seq,'ACC_ID','Prot_Seq')]
        FakeSite += [res]
phosphoBack = pd.DataFrame({'Accession' : Accessions,
                              'Sequence' : Sequences,
                              'FakeSite' : FakeSite})

#find the RK 7-11 top hit patterns for both background and actual sites

phosphosite['RK_C7-11'] =   phosphosite.apply(
            lambda row: AAwindow(row['Prot_Seq'], row['Site'],Window = (7,12), aa = 'RK'),
                           axis=1)
phosphoBack['RK_C7-11'] =   phosphoBack.apply(
            lambda row: AAwindow(row['Sequence'], row['FakeSite'],Window = (7,12), aa = 'RK'),
                           axis=1)

#what is the distribution relative to background (ie how do we get the histogram in Figure 4E)

# work out background distributions.
n = 10**5
backg = []
for i in range(n):
    wherei = np.random.choice(phosphoBack.index,1817,replace = False)
    backg += [wherei]
backg_np = np.array(backg)
tempBack = phosphoBack['RK_C7-11'].values[backg_np]
backHist = (np.count_nonzero(tempBack,axis =1))/float(1817)
phosphoRK7_11 = float(sum(phospho['RK_C7-12']>0))/1817
ax = plt.hist(backHist,bins =30,density = True)
ax = plt.vlines(x = phosphoRK7_11, ymin = 0, ymax =35,  color = 'red')
print float(sum(backHist > phosphoRK7_11))/10**5 # for writing actual number
plt.savefig('percentRK_morethan0.pdf')
plt.clf()


#find aa patterns for 1, 3, and 7 aa window 


aas = ['A','C','D','D','N','E','DE','E','Q','F','W','Y','G','H',
      'I','V','L','K','KR','R','M','P','S','ST','T']

cols1 = []
for i in range(1,30):
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i)
        cols1 += [Col]
        phosphosite[Col] = phosphosite.apply(
            lambda row: AAwindow(row['Prot_Seq'], row['Site'],Window = (i,i+1), aa = k),
                               axis=1)
        phosphoBack[Col] = phosphoBack.apply(
            lambda row: AAwindow(row['Sequence'], row['FakeSite'],Window = (i,i+1), aa = k),
                               axis=1)
cols3 = []
for i in range(1,28):
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i+2)
        cols3 += [Col]
        phosphosite[Col] = phosphosite.apply(
            lambda row: AAwindow(row['Prot_Seq'], row['Site'],Window = (i,i+3), aa = k),
                               axis=1)
        phosphoBack[Col] = phosphoBack.apply(
            lambda row: AAwindow(row['Sequence'], row['FakeSite'],Window = (i,i+3), aa = k),
                               axis=1)
cols7 = []
for i in range(1,24):
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i+6)
        cols7 += [Col]
        phosphosite[Col] = phosphosite.apply(
            lambda row: AAwindow(row['Prot_Seq'], row['Site'],Window = (i,i+7), aa = k),
                               axis=1)
        phosphoBack[Col] = phosphoBack.apply(
            lambda row: AAwindow(row['Sequence'], row['FakeSite'],Window = (i,i+7), aa = k),
                               axis=1)





out = {}
for j,k in zip(['1','3','7'],[cols1,cols3,cols7]):
    out[j] = {'Col' : k}
    out[j]['back_nonzero'] = np.zeros(len(k))
    out[j]['actual_nonzero'] = np.zeros(len(k))
    out[j]['nonzero_pval'] = np.zeros(len(k))
    for ii,i in enumerate(k):
        tempBack = phosphoBack[i].values[backg_np]
        out[j]['back_nonzero'][ii] = float(np.count_nonzero(tempBack))/(1817*n)
        out[j]['actual_nonzero'][ii] = float(np.count_nonzero(phospho[i].values))/(1817)
        if (out[j]['actual_nonzero'][ii] > out[j]['back_nonzero'][ii]):
            out[j]['nonzero_pval'][ii] = 2 * float(sum(
                ((np.count_nonzero(tempBack,axis =1))/float(1817)) > 
                                                       out[j]['actual_nonzero'][ii])
                                                      )/n
        else:
            out[j]['nonzero_pval'][ii] = 2 * float(sum(
                ((np.count_nonzero(tempBack,axis =1))/float(1817)) < 
                                                           out[j]['actual_nonzero'][ii])
                                                      )/n


#split up and do multiple testing correction
out1 = pd.DataFrame(out['1'])
out3 = pd.DataFrame(out['3'])
out7 = pd.DataFrame(out['7'])

out1.loc[out1.nonzero_pval == 0, 'nonzero_pval'] = 10**-5 # limits of analysis, so log =/= inf
out1['nonzero_padj'] = multipletests(out1.nonzero_pval)[1]
out1['pNonzero'] = -np.log10(out1.nonzero_padj)

out3.loc[out3.nonzero_pval == 0, 'nonzero_pval'] = 10**-5 # limits of analysis, so log =/= inf
out3['nonzero_padj'] = multipletests(out3.nonzero_pval)[1]
out3['pNonzero'] = -np.log10(out3.nonzero_padj)

out7.loc[out7.nonzero_pval == 0, 'nonzero_pval'] = 10**-5 # limits of analysis, so log =/= inf
out7['nonzero_padj'] = multipletests(out7.nonzero_pval)[1]
out7['pNonzero'] = -np.log10(out7.nonzero_padj)

out1.to_csv('out1.csv',index = False)
out3.to_csv('out3.csv',index = False)
out7.to_csv('out7.csv',index = False)

#get positions (avg)
out1['position'] = 0
out3['position'] = 0
out7['position'] = 0
​
​
​
aas = ['A','C','D','N','E','DE','Q','F','Y','W','G','H',
      'I','V','L','K','KR','R','M','P','S','ST','T']
cols1 = []
for i in range(4,30):
    cols1i = []
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i)
        cols1 += [Col]
        cols1i += [Col]
    out1.loc[out1.Col.isin(cols1i),'position'] = i
cols3 = []
for i in range(4,28):
    cols3i = []
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i+2)
        cols3 += [Col]
        cols3i += [Col]
    out3.loc[out3.Col.isin(cols3i),'position'] = i + 1
​
​
for i in range(4,24):
    cols7i = []
    for k in aas:
        Col = k + '_'+'C'+str(i)+'-'+str(i+6)
        cols7 += [Col]
        cols7i += [Col]
    out7.loc[out7.Col.isin(cols7i),'position'] = i + 3
    
#Plot aa
currdata = out1.reset_index(drop = True)
colors = ['#6699CC','#CC6677','#AA4499',
          '#999933',
          '#DDCC77','#117733',
         '#44AA99','#000000','#111111','#202020']
AAs = {'RK' : currdata.Col.str.contains(r'^KR_C')}
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^K_C')
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^R_C')
AAs['DE'] = currdata.Col.str.contains(r'^DE_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^D_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^E_C')
AAs['NQ'] = currdata.Col.str.contains(r'^N_C')
AAs['NQ'] = AAs['NQ'] | currdata.Col.str.contains(r'^Q_C')
AAs['IL'] = currdata.Col.str.contains(r'^I_C')
AAs['IL'] = AAs['IL'] | currdata.Col.str.contains(r'^L_C')
AAs['FWY'] = currdata.Col.str.contains(r'^Y_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^F_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^W_C')
AAs['ST'] = currdata.Col.str.contains(r'^ST_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^S_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^T_C')
AAs['VA'] = currdata.Col.str.contains(r'^V_C')
AAs['VA'] = AAs['VA'] | currdata.Col.str.contains(r'^A_C')
AAs['MGC'] = currdata.Col.str.contains(r'^C_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^M_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^G_C')
AAs['P'] = currdata.Col.str.contains(r'^P_C')
AAs['H'] = currdata.Col.str.contains(r'^H_C')
​
​
nothits1 = (currdata.nonzero_padj >= 0.05) | (np.absolute(currdata.nonzero_diff) < 0.05)
​#note- commented-out lines below are due to lack of data meeting hit criteria
ax = plt.gca()
currdata.loc[nothits1,:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = '#DDDDDD',
                         ax = ax)
# currdata.loc[AAs['RK']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[0],
#                          ax = ax)
currdata.loc[AAs['DE']& (~(nothits1)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[1],
                         ax = ax)
# currdata.loc[AAs['NQ']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[2],
#                          ax = ax)
# currdata.loc[AAs['IL']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[3],
#                          ax = ax)
# currdata.loc[AAs['FWY']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[4],
#                          ax = ax)
currdata.loc[AAs['ST']& (~(nothits1)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[5],
                         ax = ax)
# currdata.loc[AAs['VA']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[6],
#                          ax = ax)
# currdata.loc[AAs['P']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[7],
#                          ax = ax)
# currdata.loc[AAs['MGC']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[8],
#                          ax = ax)
# currdata.loc[AAs['H']& (~(nothits1)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[9],
#                          ax = ax)
plt.title('1 Amino Acid enrichments')
plt.xlabel('Position relative to site of glycosylation')
plt.ylabel('Difference % >0 /n(Glycosites - Random Ser/Thr)')
plt.legend(['nonsig','DE','ST'])
plt.savefig('./position_distance_1aa_nonzero.pdf')
plt.clf()


#plot for 3 aa windows
currdata = out3.reset_index(drop = True)
colors = ['#6699CC','#CC6677','#AA4499',
          '#999933',
          '#DDCC77','#117733',
         '#44AA99','#000000','#111111','#202020']
AAs = {'RK' : currdata.Col.str.contains(r'^KR_C')}
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^K_C')
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^R_C')
AAs['DE'] = currdata.Col.str.contains(r'^DE_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^D_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^E_C')
AAs['NQ'] = currdata.Col.str.contains(r'^N_C')
AAs['NQ'] = AAs['NQ'] | currdata.Col.str.contains(r'^Q_C')
AAs['IL'] = currdata.Col.str.contains(r'^I_C')
AAs['IL'] = AAs['IL'] | currdata.Col.str.contains(r'^L_C')
AAs['FWY'] = currdata.Col.str.contains(r'^Y_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^F_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^W_C')
AAs['ST'] = currdata.Col.str.contains(r'^ST_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^S_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^T_C')
AAs['VA'] = currdata.Col.str.contains(r'^V_C')
AAs['VA'] = AAs['VA'] | currdata.Col.str.contains(r'^A_C')
AAs['MGC'] = currdata.Col.str.contains(r'^C_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^M_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^G_C')
AAs['P'] = currdata.Col.str.contains(r'^P_C')
AAs['H'] = currdata.Col.str.contains(r'^H_C')
nothits3 = (currdata.nonzero_padj >= 0.05) | (np.absolute(currdata.nonzero_diff) < 0.05)
​
ax = plt.gca()
plotLines(currdata, nothits3,'position','nonzero_diff','colWidth',1,'#DDDDDD',
                         ax = ax)
currdata.loc[nothits3,:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = '#DDDDDD',
                         ax = ax)
​
plotLines(currdata, (AAs['P']& (~(nothits3))),'position',
          'nonzero_diff','colWidth',1,colors[7],
                         ax = ax)
currdata.loc[AAs['P']& (~(nothits3)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[7],
                         ax = ax)
plotLines(currdata, (AAs['DE']& (~(nothits3))),'position',
          'nonzero_diff','colWidth',1,colors[1],
                         ax = ax)
currdata.loc[AAs['DE']& (~(nothits3)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[1],
                         ax = ax)
​
plotLines(currdata, (AAs['IL']& (~(nothits3))),'position',
          'nonzero_diff','colWidth',1,colors[3],
                         ax = ax)
currdata.loc[AAs['IL']& (~(nothits3)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[3],
                         ax = ax)
plotLines(currdata, (AAs['ST']& (~(nothits3))),'position',
          'nonzero_diff','colWidth',1,colors[5],
                         ax = ax)
currdata.loc[AAs['ST']& (~(nothits3)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[5],
                         ax = ax)
# plotLines(currdata, (AAs['VA']& (~(nothits3))),'position',
#           'nonzero_diff','colWidth',3,colors[4],
#                          ax = ax)
# currdata.loc[AAs['VA']& (~(nothits3)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[4],
#                          ax = ax)
# plotLines(currdata, (AAs['RK']& (~(nothits3))),'position',
#           'nonzero_diff','colWidth',3,colors[0],
#                          ax = ax)
# currdata.loc[AAs['RK']& (~(nothits3)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[0],
#                          ax = ax)
# plotLines(currdata, (AAs['NQ']& (~(nothits3))),'position',
#           'nonzero_diff','colWidth',3,colors[2],
#                          ax = ax)
# currdata.loc[AAs['NQ']& (~(nothits3)),:].plot(
#     x='position',y='nonzero_diff',kind = 'scatter',color = colors[2],
#                          ax = ax)
​
plt.title('3 Amino Acid enrichments')
plt.xlabel('Position relative to site of glycosylation')
plt.ylabel('Difference % >0 (Glycosites - Random Ser/Thr)')
plt.savefig('3aa_withbars.pdf')
plt.clf()

currdata = out7.reset_index(drop = True)
colors = ['#6699CC','#CC6677','#AA4499',
          '#999933',
          '#DDCC77','#117733',
         '#44AA99','#000000','#111111','#202020']
AAs = {'RK' : currdata.Col.str.contains(r'^KR_C')}
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^K_C')
AAs['RK'] = AAs['RK'] | currdata.Col.str.contains(r'^R_C')
AAs['DE'] = currdata.Col.str.contains(r'^DE_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^D_C')
AAs['DE'] = AAs['DE'] | currdata.Col.str.contains(r'^E_C')
AAs['NQ'] = currdata.Col.str.contains(r'^N_C')
AAs['NQ'] = AAs['NQ'] | currdata.Col.str.contains(r'^Q_C')
AAs['IL'] = currdata.Col.str.contains(r'^I_C')
AAs['IL'] = AAs['IL'] | currdata.Col.str.contains(r'^L_C')
AAs['FWY'] = currdata.Col.str.contains(r'^Y_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^F_C')
AAs['FWY'] = AAs['FWY'] | currdata.Col.str.contains(r'^W_C')
AAs['ST'] = currdata.Col.str.contains(r'^ST_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^S_C')
AAs['ST'] = AAs['ST'] | currdata.Col.str.contains(r'^T_C')
AAs['VA'] = currdata.Col.str.contains(r'^V_C')
AAs['VA'] = AAs['VA'] | currdata.Col.str.contains(r'^A_C')
AAs['MGC'] = currdata.Col.str.contains(r'^C_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^M_C')
AAs['MGC'] = AAs['MGC'] | currdata.Col.str.contains(r'^G_C')
AAs['P'] = currdata.Col.str.contains(r'^P_C')
AAs['H'] = currdata.Col.str.contains(r'^H_C')
nothits7 = (currdata.nonzero_padj >= 0.05) | (np.absolute(currdata.nonzero_diff) < 0.05)
​
ax = plt.gca()
plotLines(currdata, nothits7,'position','nonzero_diff','colWidth',3,'#DDDDDD',
                         ax = ax)
currdata.loc[nothits7,:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = '#DDDDDD',
                         ax = ax)
​
plotLines(currdata, (AAs['P']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[7],
                         ax = ax)
currdata.loc[AAs['P']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[7],
                         ax = ax)
plotLines(currdata, (AAs['DE']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[1],
                         ax = ax)
currdata.loc[AAs['DE']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[1],
                         ax = ax)
​
plotLines(currdata, (AAs['IL']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[3],
                         ax = ax)
currdata.loc[AAs['IL']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[3],
                         ax = ax)
plotLines(currdata, (AAs['ST']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[5],
                         ax = ax)
currdata.loc[AAs['ST']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[5],
                         ax = ax)
plotLines(currdata, (AAs['VA']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[4],
                         ax = ax)
currdata.loc[AAs['VA']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[4],
                         ax = ax)
plotLines(currdata, (AAs['RK']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[0],
                         ax = ax)
currdata.loc[AAs['RK']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[0],
                         ax = ax)
plotLines(currdata, (AAs['NQ']& (~(nothits7))),'position',
          'nonzero_diff','colWidth',3,colors[2],
                         ax = ax)
currdata.loc[AAs['NQ']& (~(nothits7)),:].plot(
    x='position',y='nonzero_diff',kind = 'scatter',color = colors[2],
                         ax = ax)
​
plt.title('7 Amino Acid enrichments')
plt.xlabel('Position relative to site of glycosylation')
plt.ylabel('Difference % >0 (Glycosites - Random Ser/Thr)')
plt.savefig('7aa_withbars.pdf')
plt.clf()

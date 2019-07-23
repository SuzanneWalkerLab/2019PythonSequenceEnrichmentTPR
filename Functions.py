#!/usr/bin/env python2

#a document of all new functions used in python code for Joiner, Levine et al.
import requests
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

#for ReadIn.py

def GetSeqUniprot(Accession):
    if '; ' in Accession:
        Accession = Accession.split('; ')[0]
    result = requests.get(BASE + KB_ENDPOINT + Accession + '.fasta')
    if result.ok:
        x = result.text.split('\n')
        x = ''.join(x[1:])
    else:
        print("Something went wrong: ", result.status_code)
    return x

def whereLower(inputString):
    for i in range(0,len(inputString)):
        if(inputString[i].islower() & ((inputString[i] == 't') | (inputString[i] == 's'))):
            output = i
            break
    return output

def pepFind(Seq):
    '''place Sequence in first position, Peptide in second'''
    lfind = Seq[0].find(Seq[1].upper())
    rfind = Seq[0].rfind(Seq[1].upper())
    if rfind!=lfind:
        x = 0
    else:
        x = lfind
    return x

def getSite(Seq, Nterm, Cterm):
    '''Find location of glycosylation and return N-term amino acids before glycosite and C-term amino acids after,
    with "-" for end of protein'''
    if (Seq[1]-Nterm)<0:
        start = '-'*abs(Seq[1]-Nterm) + Seq[0][Seq[1]-(Nterm-abs(Seq[1]-Nterm)):Seq[1]]
    else:
        start = Seq[0][Seq[1]-(Nterm):Seq[1]]
    site = Seq[0][Seq[1]].lower()
    if (Cterm+Seq[1])>len(Seq[0]):
        over = Seq[1]+Cterm-len(Seq[0])
        end =  Seq[0][Seq[1]+1:Seq[1]+Cterm-over] + '-'*over
    else:
        end = Seq[0][Seq[1]+1:Seq[1]+Cterm] 
    return start + site + end

def SiteDef(data):
    #from having Site of peptide, introduce site
    data['Site'] = 0
    data['multiSite'] = False
    data = data.astype('object')
    for i in data.index:
        z = data.loc[i,'Site_Seq']
        if isinstance(data.loc[i,'SitePep'],list):
            data.at[i,'Site'] = [x + z for x in data.loc[i,'SitePep']]
            data.loc[i,'multiSite'] = True
        else:
            data.loc[i,'Site'] =  data.loc[i,'SitePep'] + z
    return data

def findMultiSite(data):
    #finds where peptide does NOT uniquely map, returns indices
    multiSite = []
    for i in data.index:
        if isinstance(data.loc[i,'Site'],list):
            multiSite = multiSite + [i]
    return multiSite

def removeSiteDuplicates(data):
    #there are, in some instances, multiple peptides for a given site. This function adds the numbers in "WT"
    #and "Mutant" to the first occurence and drops the other rows for the same site
    for Site in data.Site.unique():   
        if len(data.Accession[data.Site == Site].unique()) < sum(data.Site == Site):
            for Accession in data.Accession[data.Site == Site].unique():
                whereSite = (data.Site == Site) & (data.Accession == Accession)
                whichSite = [i for i, x in enumerate(whereSite) if x]
                data.loc[whichSite[0],'WT'] = data.WT[whereSite].sum()
                data.loc[whichSite[0],'Mutant'] = data.Mutant[whereSite].sum()
                data = data.drop_duplicates(subset = ['Accession','Site'], keep = 'first').reset_index(drop = True)
    data = data.drop(columns = ['Peptide','Site_Seq','SitePep'])
    data['Site_Seq'] = [x[1]['Sequence'][x[1]['Site']-3:x[1]['Site']+50] for x in data.iterrows()]
    return data
def mergeData(data1,data2,repInfo = True):
    #takes in two data sets and collapses to one
    if repInfo:
        data1['WT1'] = data1.WT
        data1['Mutant1'] = data1.Mutant
        data1['WT2'] = 0
        data1['Mutant2'] = 0
        data2['WT1'] = 0
        data2['Mutant1'] = 0
        data2['WT2'] = data2.WT
        data2['Mutant2'] = data2.Mutant
    for i in data1.index:
        #find the same protein in replicate1  
        if np.any(data2.Site[data2.Accession == data1.Accession[i]] == data1.Site[i]):
            data2loc = (data2.Accession == data1.Accession[i]) & (data2.Site == data1.Site[i])
            data2.loc[data2loc,'WT'] += data1.WT[i]
            data2.loc[data2loc,'Mutant'] += data1.Mutant[i]
            if repInfo:
                data2.loc[data2loc,'WT1'] = data1.WT1[i]
                data2.loc[data2loc,'Mutant1'] = data1.Mutant1[i]
        else:
            data2 = data2.append(data1.loc[i,:].copy(),ignore_index = True).reset_index(drop = True)
    return data2


#functions to export FASTA for WebLogos

def getSite(Seq, Nterm, Cterm):
    '''Find location of glycosylation and return N-term amino acids before glycosite and C-term amino acids after,
    with "-" for end of protein'''
    if (Seq[1]-Nterm)<0:
        start = '-'*abs(Seq[1]-Nterm) + Seq[0][Seq[1]-(Nterm-abs(Seq[1]-Nterm)):Seq[1]]
    else:
        start = Seq[0][Seq[1]-(Nterm):Seq[1]]
    site = Seq[0][Seq[1]].lower()
    if (Cterm+Seq[1])>len(Seq[0]):
        over = Seq[1]+Cterm-len(Seq[0])
        end =  Seq[0][Seq[1]+1:Seq[1]+Cterm-over] + '-'*over
    else:
        end = Seq[0][Seq[1]+1:Seq[1]+Cterm] 
    return start + site + end

def toFASTA(SeqDF, identifierCol, SeqCol):
    '''Takes in a pair of pandas columns, an identifier and a sequence, and outputs as fasta string'''
    x = ''
    for i in range(0,SeqDF.shape[0]):
        x = x + '>'+ str(i) + SeqDF[identifierCol][i]  + '\n' + SeqDF[SeqCol][i] + '\n'
    return x

#analysis of enrichment between OGT variants in HeLa extracts

def AAwindow(Sequence,Site,Window,aa):
    #looks at sequence around site based on Window (a tuple of starting and ending
    #range, which and checks how many  of aa listed in aa are in the window. Note: Window
    #indexes site at 0
    howMany = 0
    start = Window[0]+Site-1#-1 to deal with 0-indexing
    if start <0:
        start = 0
    for i in Sequence[start:Window[1]+Site-1]:#note non-inclusive slicing of second number in python
        if i in aa:
            howMany += 1
    return howMany

def mannUpval(cols, data, comparisons, compDesc, labels, mean = True, median = False,
              fold = True, diff = True, storeU = False):
    #systematize this...
    #data is a list of dataframes, comparisons is a list saying
    #indices of dataframes to compare in data, and compDesc gives the string names to use for the comparison
    output = {'Cols':cols}
    for ik, k in enumerate(labels):        
        if mean:
            means  = []
        if median:
            medians = []
        for j in cols:
            if mean:
                means  += [data[ik][j].mean()]
            if median:
                medians += [data[ik][j].median()]
        if mean:
            output[k + '_mean' ]  = means
        if median:
            output[k + '_median' ] = medians
    for ik, k in enumerate(compDesc):
        pvals = []
        if storeU:
            Uvals = []
        if fold:
            folds = []
        if diff:
            diffs = []
        for i in cols:
            if (sum(data[comparisons[ik][0]][i]>0) + 
                sum(data[comparisons[ik][1]][i]>0)) > 20:
                try:
                    x,y = mannwhitneyu(data[comparisons[ik][0]][i],
                                       data[comparisons[ik][1]][i],alternative = 'two-sided')
                except:
                    x,y = np.nan,np.nan
            else:
                x,y = np.nan,np.nan
            pvals += [y]
            if storeU:
                Uvals += [x]
            if fold:
                try:
                    folds += [float(data[comparisons[ik][0]][i].mean())/float(
                        data[comparisons[ik][1]][i].mean())]
                except:
                    folds += [np.nan]
            if diff:
                diffs += [float(data[comparisons[ik][0]][i].mean()) - float(
                    data[comparisons[ik][1]][i].mean())]
        output[k+'_pval'] = pvals
        output['p' + k] = -np.log10(pvals)
        if storeU:
            output[k + '_U'] = Uvals
        if fold:
            output[k + '_fold'] = folds
            output[k + '_log2fold'] = np.log2(output[k + '_fold'])
        if diff:
            output[k + '_diff'] = diffs
    return pd.DataFrame(output)  

#for phosphosite analysis

def findST(Sequence): #for generating background distribution
    #finds all S/T occurrences
    ret = []
    for ix,x in enumerate(Sequence):
        if x in ['S','T']:
            ret += [ix+1]
    return ret
# for the plots of windows in phosphosite data
def plotLines(data, rows, columnx, columny, columnWidth, padding, color, ax):
    data = data.reset_index(drop = True)
    dataToUse1 = data.loc[rows,:].reset_index(drop = True)
    dataToUse2 = dataToUse1.append(dataToUse1,ignore_index = True)
    for ii, i in enumerate(dataToUse1[columnx]):
        x = dataToUse2.loc[dataToUse2.Col == dataToUse2.Col[ii],:].reset_index(drop = True)
        x.loc[0,columnx] = x.loc[0,columnx]-padding
        x.loc[1,columnx] = x.loc[1,columnx]+padding
        x.plot(x=columnx,y=columny,kind = 'line',color = color,
                         ax = ax, lw = 0.8, 
               alpha = 0.2, legend = False)
​

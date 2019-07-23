#!/usr/bin/env python2

#reads in data from .csv's derived from raw data,
#the column headings are "Description" (protein description,
# "Accesion" (UniProt Accession),
# "RowLabels" (Tryptic peptide sequence with lowercase site of modification; note we change this to "Peptide" early on)
# "WildType" (number of times a peptide was detected in this replicate for WT-treated Extracts)
# and "Mutant" (number of times a peptide was detected in this replicate for D2A-treated Extracts)

import numpy as np
import pandas as pd
import requests
from Functions.py import GetSeqUniprot
from Functions.py import whereLower
from Functions.py import pepFind
from Functions.py import removeSiteDuplicates
from Functions.py import mergeData

replicate1 = pd.read_csv('DataReplicate1.csv')
replicate2 = pd.read_csv('DataReplicate2.csv')

BASE = "http://www.uniprot.org"
KB_ENDPOINT = "/uniprot/"


replicate1.columns = ['Description','Accession','Peptide','WT','Mutant']
replicate2.columns = ['Description','Accession','Peptide','WT','Mutant']
replicate1['Sequence'] = replicate1.Accession.apply(GetSeqUniprot)
replicate2['Sequence'] = replicate2.Accession.apply(GetSeqUniprot)
replicate1['Site_Seq'] = replicate1.Peptide.apply(whereLower)
replicate1['SitePep'] = replicate1[['Sequence','Peptide']].apply(pepFind,axis = 1)
replicate1['Site'] = replicate1.SitePep + replicate1.Site_Seq
replicate2['Site_Seq'] = replicate2.Peptide.apply(whereLower)
replicate2['SitePep'] = replicate2[['Sequence','Peptide']].apply(pepFind,axis = 1)
replicate2['Site'] = replicate2.SitePep + replicate2.Site_Seq
#this block takes care of not-uniquely-mapping sites
replicate1 = SiteDef(replicate1)
replicate2 = SiteDef(replicate2)
replicate1 = replicate1.drop(index = findMultiSite(replicate1)).reset_index(drop = True)
replicate2 = replicate2.drop(index = findMultiSite(replicate2)).reset_index(drop = True)
replicate1 = replicate1.drop(columns = 'multiSite')#removes marker column
replicate2 = replicate2.drop(columns = 'multiSite')#removes marker column
#merge replicated sites (ie ones with multiple different peptides corresponding to same site)
replicate1 = removeSiteDuplicates(replicate1).sort_values(by = 'Accession').reset_index(drop = True)
replicate2 = removeSiteDuplicates(replicate2).sort_values(by = 'Accession').reset_index(drop = True)  
data = mergeData(replicate1,replicate2)
data.to_csv('data.csv',index=False)
replicate2.to_csv('replicate2.csv',index = False)
replicate1.to_csv('replicate1.csv',index = False)



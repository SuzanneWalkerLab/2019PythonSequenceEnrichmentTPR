#!/usr/bin/env python2

#takes data in various forms and outputs sequence directly around site as FASTA in alignment format for direct input to WebLogo3

import pandas as pd
from Functions.py import getSite
from Functions.py import toFASTA

data = pd.read_csv('data.csv', header = 0)
data['Glyco15'] = data[['Sequence','Site']].apply(getSite,Nterm=0,Cterm=15,axis = 1)
WT = data.WT > 0
WTonly = WT & (data.Mutant == 0)
D2A = data.Mutant > 0
D2Aonly = D2A & (data.WT == 0)
WT_all_fasta  = toFASTA(data.loc[WT, :].reset_index(),'Accession','Glyco15')
WT_specific_fasta  = toFASTA(data.loc[WTonly, :].reset_index(),'Accession','Glyco15')
D2A_all_fasta  = toFASTA(data.loc[D2A, :].reset_index(),'Accession','Glyco15')
D2A_specific_fasta  = toFASTA(data.loc[D2Aonly, :].reset_index(),'Accession','Glyco15')
#output fasta for each
with open("WT_all.fasta", "w") as text_file:
    text_file.write(WT_all_fasta)
with open("D2A_all.fasta", "w") as text_file:
    text_file.write(D2A_all_fasta)
with open("WT_specific.fasta", "w") as text_file:
    text_file.write(WT_specific_fasta)
with open("D2A_specific.fasta", "w") as text_file:
    text_file.write(D2A_specific_fasta)

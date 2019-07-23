#!/usr/bin/env python2
#figure out site statistics on protein and glycosite level (related to Figure 3B, 3C, S8)

import pandas as pd
import numpy as np

data = pd.read_csv('data.csv')

#define the categories

WTonly = (data.WT > 0) & (data.Mutant == 0)
D2Aonly = (data.Mutant > 0) & (data.WT == 0)
Both = (data.WT > 0) & (data.Mutant > 0)

#site data (related to figure 3C)
print "WT only site number", sum(WTonly)
print "D2A only site number", sum(D2Aonly)
print "WT and D2A site number", sum(Both)

#protein data
WTproteins = data.Accession[data.WT > 0].unique()
D2Aproteins = data.Accession[data.Mutant > 0].unique()
BothLogicalWT = np.isin(WTproteins,D2Aproteins)
WTonlyProtLog = np.invert(BothLogicalWT)
D2AonlyProtLog = np.invert(np.isin(D2Aproteins,WTproteins))
print "WT only protein number", sum(WTonlyProtLog)
print "D2A only protein number", sum(D2AonlyProtLog)
print "WT and D2A protein number", sum(BothLogicalWT)

#data for sites on "Both" proteins
bothData = data.loc[data.Accession.isin(WTproteins[BothLogicalWT]),].reset_index(drop = True)
bothWTonly = (bothData.WT > 0) & (bothData.Mutant == 0)
bothD2Aonly = (bothData.Mutant > 0) & (bothData.WT == 0)
bothBoth = (bothData.WT > 0) & (bothData.Mutant > 0)
print "WT only site number (proteins glycosylated by both)", sum(bothWTonly)
print "D2A only site number (proteins glycosylated by both)", sum(bothD2Aonly)
print "WT and D2A site number (proteins glycosylated by both)", sum(bothBoth)



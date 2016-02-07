# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 22:33:32 2016

@author: Hanjo
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def fig_phosints():
    """
    Visualize the non-phospgorylated vs. phosphorylated peptide ratios.
    """
    pass

sheetnames = ["log10HL", "log10HM", "log10ML"]
for sheet in sheetnames:
    excel = pd.read_excel("D:\Sven\Dropbox\VX data for Sven\SVEN\data\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide_overview.xlsx", sheetname=sheet)
    break


f, ax = plt.subplots(1, figsize=(11.69, 8.27))


excel.columns


plt.hist(excel[excel["direction"]=="Down"])

plt.hist(excel[excel["direction"]=="Up"])
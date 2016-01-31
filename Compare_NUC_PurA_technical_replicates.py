# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 16:14:57 2015
 
@author: sven
 
 
 
Project descript goes here:
 
0) remove nans
0.1) missing data imputation
0.2)
 
1) Normalize phosphopeptides by normal peptides
2) compare the ratios of proteins and peptides
2.1) does the phosphopeptide change come from a change in proteine xpression?
2.2) or does it come by a change in phosphorylation
3)
 
 
Get everything ready
 
cd /home/sven/software/sklearn/scikit-learn/
import sklearn
cd /home/sven/workspace/git_workspace/projects/cytokinesis_2015/KNN/
#do this by hand?
%run imputation.py
"""
 
# global variables
multiple_sites = 0
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import re
import os
import brewer2mpl
import HTSeq
from scipy.stats import norm
import scipy.stats as stat
import matplotlib as mpl
import itertools
from cyto_utils import *
 
 
#    TODO
#==============================================================================
# - add phosphosite annotation
# - add phosphosite known information (also from PSP)
# -
#==============================================================================
path = "/home/sven/Dropbox/VX data for Sven/SVEN/"
#path = "E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\"
path = "/home/sven/Dropbox/VX data for Sven/SVEN/"
path = "E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\"
path = "E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\data\\"
path = "/home/sven/Dropbox/VX data for Sven/SVEN/"
path = "C:\\Users\\Han\\Dropbox\\VX data for Sven\\SVEN\\data\\"
path = "D:\\Sven\\Dropbox\\VX data for Sven\\SVEN\\"
path_input = "D:\\Sven\\Dropbox\\VX data for Sven\\technical replicates\\"
fasta_file = "C:\\Users\\Han\\Dropbox\\VX data for Sven\\sprot_2014_08_2014_11.fasta"
fasta_file = "D:\\Sven\\Dropbox\\VX data for Sven\\sprot_2014_08_2014_11.fasta"
 

files = []
for root, dirnames, filenames in os.walk(path_input):
    files.extend(glob.glob(root + "/*.xlsx"))
 
 
kinomes = []
for root, dirnames, filenames in os.walk(path):
    kinomes.extend(glob.glob(root + "/*full*.tsv"))
    
peptide_files = sorted([i for i in files if ("overview" not in i)])
#kinome_files = sorted(kinomes)
 

print "Los gehts...."
 

#files = glob.glob(path+"*.xlsx")
#files = [i for i in files if "overview" not in i]
 
#files = glob.glob(path+"*.xlsx")
#files = [i for i in files if "overview" not in i]
 
#peptide_files = sorted([i for i in files if "peptide" in i])
#protein_files = sorted([i for i in files if "protein" in i])
 

 
FASTA_dic = get_fasta_dic(fasta_file)
phosphositeDB = read_phosphosite(path+"\\Phosphosite_data\\Phosphorylation_site_dataset")
 

reports = []
 
ratios_default = ["Heavy/Light"]
ratios_log = ['log10HL']
ratios_norm = ['norm_log10HL']
alpha = 0.05
columns = [u'# Proteins', 'Sequence', u'# Protein Groups', u'Protein Group Accessions',
           u'Modifications', u'MH+ [Da]',
           u'phosphoRS Site Probabilities', u'Area',
           u'Heavy/Light Count', u'Heavy/Light Variability [%]',
           u'Medium/Light Count',
           

           u'Medium/Light Variability [%]', u'q-Value', u'PEP',
           u'# Missed Cleavages',
           u'Medium/Light', u'Heavy/Light']
           #"Heavy/Medium", "isPhospho"]

dfs_processed_phos = []
dfs_processed_reg = []
dfs_all = []
#%%
for peptides in peptide_files:    
    report = ""
    filename = os.path.basename(peptides)
    print filename
    print "Get  peptides tables....({})".format(filename)
    peptides_df = pd.ExcelFile(peptides).parse("Sheet1")
    peptides_df = peptides_df[columns]
    peptides_df["isPhospho"] = pd.notnull(peptides_df["phosphoRS Site Probabilities"])
 
    print "Filter phosphopeptides..."
    phospho_peptides = peptides_df[peptides_df["isPhospho"]==True].copy()
    regular_peptides = peptides_df[peptides_df["isPhospho"]==False].copy()
 
     
    # annotate sequence_window
    print "Annotate sequence windows..."
    add_seq_window(phospho_peptides, FASTA_dic)
    #transform the ratios to log10
    print "Compute logarithms..."
    print "\tPhosphopeptides..."
    computelog(phospho_peptides, ratio_columns=["Heavy/Light"])
    print "\tRegular Peptides..."
    computelog(regular_peptides, ratio_columns=["Heavy/Light"])
    
    #add phosphosite plus site information
    print "Add PhosphoSitePlus data"
    #add_known_pssite(phospho_peptides, phosphositeDB, FASTA_dic) 
    
    #delete non-existing columns
    del_columns = ['Medium/Light', 'Heavy/Medium', u'log10ML', u'log10HM']
    for del_i, df_i in zip(del_columns, [peptides_df, phospho_peptides, regular_peptides]):
        try:
            del df_i[del_i]
        except:
            pass
        
    dfs_processed_phos.append(phospho_peptides)
    dfs_processed_reg.append(phospho_peptides)
    
    #==========================================================================
    #     Analyze and store as single Excel file
    #==========================================================================
    dfs = []
    for ratio_i in ratios_log:
        print ratio_i
        dfs.append(analyze_ratio(phospho_peptides, regular_peptides,
                                 ratio_i, alpha, peptides[:-5]))
 
    writer = pd.ExcelWriter(peptides[:-5] + '____technical_overview.xlsx',
                            engine='xlsxwriter')
    for i, dataframe_i in enumerate(dfs):
        grouped = dataframe_i.groupby(["direction"])
        for j, dataframe_j in grouped:
            dataframe_j.to_excel(writer, 
                                 sheet_name="{}_{}".format(ratios_log[i], j))
 
        dataframe_i.to_excel(writer, sheet_name=ratios_log[i])
    writer.save()
    #get some stats
    frequency_tab = get_phosphorylated_stats(phospho_peptides)
    print "Results written to: {}".format(peptides[:-5] + '____technical_overview.xlsx')
    dfs_all.append(dfs)
 

    
dfs[0]["Area2"] = dfs[0]["Area"] / dfs[0]["Heavy/Light"]
 
plt.scatter(dfs[0]["Area"], dfs[0]["Area2"]) 
np.correlate(dfs[0]["Area"], dfs[0]["Area2"])
scipy.stats.pearsonr(dfs[0]["Area"], dfs[0]["Area2"])



#    report = """File: {}
#    Peptide matrix (all): {}
#    Protein matrix (all): {}
#
#    missing Heavy/Light features: {} / {}
#    missing Medium/Light features: {} / {}
#    missing HL or missing ML: {} / {}
#
#
#    regular peptides: {} / {}
#    phosphopeptides: {} / {}
#        S/T/Y (%): {} / {} / {}
#
#    overlap regular and phosphopeptides: {}%
#    """.format(os.path.basename(peptides), peptides_df.shape, proteins_df.shape,
#               len(HL_filter) - np.sum(HL_filter), len(HL_filter),
#               len(ML_filter) - np.sum(ML_filter), len(ML_filter),
#               len(filter_) - np.sum(filter_), len(filter_),
#               regular_peptides.shape[0], filtered_peptides_df.shape[0],
#               phospho_peptides.shape[0], filtered_peptides_df.shape[0],
#               frequency_tab["S"]*100, frequency_tab["T"]*100,
#               frequency_tab["Y"]*100, np.round(intersection_val*100, 2)
#               )
#    print report


#    break
#    explorative_violin(t)
#    explorative_scatter(phospho_peptides)

#    fobj = open("E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\exploratory_figures_tables\\{}_report.txt".format(filename), "w")
#    for line in report:
#        fobj.write(line)
#    fobj.close()
#
#    reports.append(report)
 
 
"""
peptides_df["Heavy/Light"][:5]
Out[64]:
0    0.927501
1         NaN
2    1.820355
3    0.841035
4         NaN
Name: Heavy/Light, dtype: float64
 
pd.notnull(peptides_df["Heavy/Light"][:5])
Out[65]:
0     True
1    False
2     True
3     True
4    False
Name: Heavy/Light, dtype: bool
"""

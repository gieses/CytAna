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
import numpy as np
import os
import cyto_utils as cutils

#    TODO
#==============================================================================
# - add phosphosite annotation
# - add phosphosite known information (also from PSP)
# -
#==============================================================================

#==============================================================================
# Read data and prepare data structures
#==============================================================================
path = "D:\\Sven\\Dropbox\\VX data for Sven\\SVEN\\"
fasta_file = "D:\\Sven\\Dropbox\\VX data for Sven\\sprot_2014_08_2014_11.fasta"
outpath = "D:\\Sven\\"

files = []
for root, dirnames, filenames in os.walk(path):
    files.extend(glob.glob(root + "/*.xlsx"))


kinomes = []
for root, dirnames, filenames in os.walk(path):
    kinomes.extend(glob.glob(root + "/*full*.tsv"))

peptide_files = sorted([i for i in files if (("peptide" in os.path.basename(i)) and ("overview" not in i) and ("Panther" not in i) and ("CytoPurA_CytoNoc" in i))])
protein_files = sorted([i for i in files if (("protein" in os.path.basename(i)) and ("overview" not in i) and ("Panther" not in i) and ("CytoPurA_CytoNoc" in i))])
kinome_files = sorted(kinomes)

#==============================================================================
# Start processing
#==============================================================================
FASTA_dic = cutils.get_fasta_dic(fasta_file)
phosphositeDB = cutils.read_phosphosite(path+"\\Phosphosite_data\\Phosphorylation_site_dataset")


reports = []

ratios_default = ["Heavy/Light"]
ratios_log = ['log10HL']
ratios_norm = ['norm_log10HL']
alpha = 0.05
columns = [u'# Proteins', 'Sequence', u'# Protein Groups', u'Protein Group Accessions',
           u'Modifications', u'MH+ [Da]', u'phosphoRS Isoform Probability',
           u'phosphoRS Site Probabilities', u'Area',
           u'Heavy/Light Count', u'Heavy/Light Variability [%]',
           u'Medium/Light Count',
           u'Medium/Light Variability [%]', u'q-Value', u'PEP',
           u'phosphoRS Binomial Peptide Score', u'# Missed Cleavages',
           u'Medium/Light', u'Heavy/Light']
           #"Heavy/Medium", "isPhospho"]


########### get the input data for kinomeexplorer
#only execute once
execute = False
if execute:
    for peptides, proteins in zip(peptide_files, protein_files):
        dirname = os.path.dirname(peptides)
        suffix = os.path.basename(peptides)
        peptides_df = pd.ExcelFile(peptides).parse("Sheet1")
        peptides_df["isPhospho"] = pd.notnull(peptides_df["phosphoRS Isoform Probability"])
        cutils.get_kinomexplorer_output(peptides_df[peptides_df["isPhospho"]==True],
                                 FASTA_dic,
                                 "{}\\kinome_{}.csv".format(dirname, suffix), fasta_file)
else:
    pass


#%%
for peptides, proteins, kinomes in zip(peptide_files, protein_files, kinome_files):
    kinome_df = cutils.read_kinome_explorer(kinomes)

    report = ""
    filename = os.path.basename(peptides)
    print "Get  peptides tables....({})".format(filename)
    proteins_df = pd.ExcelFile(proteins).parse("Sheet1")
    proteins_df["Heavy/Medium"] = np.divide(proteins_df["Heavy/Light"],
                                            proteins_df["Medium/Light"])

    peptides_df = pd.ExcelFile(peptides).parse("Sheet1")
    peptides_df = peptides_df[columns]
    peptides_df["Heavy/Medium"] = np.divide(peptides_df["Heavy/Light"],
                                            peptides_df["Medium/Light"])
    peptides_df["isPhospho"] = pd.notnull(peptides_df["phosphoRS Isoform Probability"])

    phospho_peptides = peptides_df[peptides_df["isPhospho"]==True].copy()
    regular_peptides = peptides_df[peptides_df["isPhospho"]==False].copy()

    # annotate sequence_window
    cutils.add_seq_window(phospho_peptides, FASTA_dic)
    #transform the ratios to log10
    cutils.computelog(phospho_peptides)
    cutils.computelog(regular_peptides)
    cutils.computelog(proteins_df)

    #add phosphosite plus site information
    cutils.add_known_pssite(phospho_peptides, phosphositeDB, FASTA_dic)

    #delete non-existing columns
    del_columns = ['Medium/Light', 'Heavy/Medium', u'log10ML', u'log10HM']
    for del_i, df_i in zip(del_columns, [peptides_df, phospho_peptides, regular_peptides]):
        try:
            del df_i[del_i]
        except:
            pass

    #==========================================================================
    #     Analyze and store as single Excel file
    #==========================================================================
    dfs_phosphor = []
    dfs_regular = []
    for ratio_i in ratios_log:
        print ratio_i
        temp_dfs = cutils.analyze_ratio(phospho_peptides, regular_peptides,
                                        ratio_i, alpha, peptides[:-5])
        dfs_phosphor.append(temp_dfs[0])
        dfs_regular.append(temp_dfs[1])

    #write the phosphopeptides to a excel sheet
    writer = pd.ExcelWriter(peptides[:-5] + '_overview.xlsx', engine='xlsxwriter')
    #write to single excel sheets in one file
    for i, dataframe_i in enumerate(dfs_phosphor):
        grouped = dataframe_i.groupby(["direction"])
        for j, dataframe_j in grouped:
            dataframe_j.to_excel(writer, sheet_name="{}_{}".format(ratios_log[i], j))

        dataframe_i.to_excel(writer, sheet_name=ratios_log[i])
    writer.save()

    #write the non-phosphorylated peptides to an extra file
    writer = pd.ExcelWriter(peptides[:-5] + '_nonhospho_overview.xlsx', engine='xlsxwriter')
    #write to single excel sheets in one file
    for i, dataframe_i in enumerate(dfs_regular):
        dataframe_i.to_excel(writer, sheet_name=ratios_log[i])
    writer.save()

    #get some stats
    frequency_tab = cutils.get_phosphorylated_stats(phospho_peptides)


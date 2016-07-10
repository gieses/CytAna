# -*- coding: utf-8 -*-
"""
Created on Fri May 27 06:32:13 2016

@author: Hanjo
"""

# global variables
multiple_sites = 0
import glob
import pandas as pd
import numpy as np
import os
import cyto_utils as cutils
import sys
import matplotlib.pyplot as plt
import re
import scipy.stats as stat
from statsmodels.sandbox.stats.multicomp import multipletests
#sys.exit()

#    TODO
#==============================================================================
# - add phosphosite annotation
# - add phosphosite known information (also from PSP)
# -
#==============================================================================

#==============================================================================
# Read data and prepare data structures
#==============================================================================
path = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\\May_2016_new cytokinesis files\\"
fasta_file = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\sprot_2014_08_2014_11.fasta"
outpath = "D:\\Sven\\OezgeNew\\"


path = "/home/sven/Dropbox/shared_folders/nurhan/VX data for Sven/SVEN/May_2016_new cytokinesis files/"
fasta_file = "/home/sven/Dropbox/shared_folders/nurhan/VX data for Sven/sprot_2014_08_2014_11.fasta"
outpath = "/home/sven/data/nurhan/oezge/"

files = []
for root, dirnames, filenames in os.walk(path):
    files.extend(glob.glob(root + "/*.xlsx"))


kinomes = []
for root, dirnames, filenames in os.walk(path):
    kinomes.extend(glob.glob(root + "/*full*.tsv"))

peptide_files = sorted([i for i in files if (("peptide" in os.path.basename(i)) and ("overview" not in i) and ("Panther" not in i))])
protein_files = sorted([i for i in files if (("protein" in os.path.basename(i)) and ("overview" not in i) and ("Panther" not in i))])
kinome_files = sorted(kinomes)

#==============================================================================
# Start processing
#==============================================================================
FASTA_dic = cutils.get_fasta_dic(fasta_file)
#phosphositeDB = cutils.read_phosphosite("D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\\Phosphosite_data\\Phosphorylation_site_dataset")
phosphositeDB = cutils.read_phosphosite("/home/sven/Dropbox/shared_folders/nurhan/VX data for Sven/SVEN/Phosphosite_data/Phosphorylation_site_dataset")


reports = []

ratios_default = ["Heavy/Light"]
ratios_log = ['log2HL']
ratios_norm = ['norm_2HL']
alpha = 0.05
columns = [u'# Proteins', 'Sequence', u'# Protein Groups', u'Protein Group Accessions',
           u'Modifications', u'MH+ [Da]', u'phosphoRS Isoform Probability',
           u'phosphoRS Site Probabilities', u'Area',
           u'Heavy/Light Count', u'Heavy/Light Variability [%]',
           u'q-Value', u'PEP',
           u'phosphoRS Binomial Peptide Score', u'# Missed Cleavages',
           u'Heavy/Light']
           #"Heavy/Medium", "isPhospho"]


########### get the input data for kinomeexplorer
#only execute once to
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



#store dataframes
order = []
dfs_norm_phos = []
dfs_norm_regu = []
thresholds = []
invert_ar = [False, True]
#%%
for peptides, proteins, kinomes, invert in zip(peptide_files, protein_files, kinome_files, invert_ar):
    kinome_df = cutils.read_kinome_explorer(kinomes)


    report = ""
    filename = os.path.basename(peptides)
    print "Get  peptides tables....({})".format(filename)
    proteins_df = pd.ExcelFile(proteins).parse("Sheet1")


    peptides_df = pd.ExcelFile(peptides).parse("Sheet1")
    peptides_df = peptides_df[columns]


    peptides_df["isPhospho"] = pd.notnull(peptides_df["phosphoRS Isoform Probability"])

    phospho_peptides = peptides_df[peptides_df["isPhospho"]==True].copy()
    regular_peptides = peptides_df[peptides_df["isPhospho"]==False].copy()

    if invert:
        print "inverting ratios...."
        phospho_peptides["Heavy/Light"] = 1 / phospho_peptides["Heavy/Light"]
        regular_peptides["Heavy/Light"] = 1 / regular_peptides["Heavy/Light"]

    # annotate sequence_window
    cutils.add_seq_window(phospho_peptides, FASTA_dic)
    #transform the ratios to log10
    cutils.computelog(phospho_peptides, ratio_columns=["Heavy/Light"])
    cutils.computelog(regular_peptides, ratio_columns=["Heavy/Light"])
    cutils.computelog(proteins_df, ratio_columns=["Heavy/Light"])

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
        thresholds.append(temp_dfs[2])


    #append data to combine
    dfs_norm_phos.append(dfs_phosphor)
    dfs_norm_regu.append(dfs_regular)
    order.append(os.path.basename(peptides))

    ##########################################################################

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



#==============================================================================
# here starts the combinin of the two files
#==============================================================================
df_phos1 = dfs_norm_phos[0][0].copy()
df_phos2 = dfs_norm_phos[1][0].copy()


df_reg1 = dfs_norm_regu[0][0].copy()
df_reg2 = dfs_norm_regu[1][0].copy()

plt.boxplot([df_phos1["norm_log2HL"],
             df_phos2["norm_log2HL"],
             df_reg1["norm_log2HL"],
             df_reg2["norm_log2HL"]])

plt.hist(df_phos1["norm_log2HL"], bins=30, normed=True, alpha=.7)
plt.hist(df_phos2["norm_log2HL"], bins=30, normed=True, alpha=.6)

#add identifier by phospho
df_phos1["phos_id"] = ["{}_{}".format(i, get_phospho(j)) for i,j in
                       zip(df_phos1["Sequence"], df_phos1["Modifications"])]

df_phos2["phos_id"] = ["{}_{}".format(i, get_phospho(j)) for i,j in
                       zip(df_phos2["Sequence"], df_phos2["Modifications"])]

#combine data sets
all_ids = np.union1d(df_phos1["phos_id"].values, df_phos2["phos_id"].values)

df_phos1 = df_phos1.set_index("phos_id")
df_phos2 = df_phos2.set_index("phos_id")
print thresholds


#==============================================================================
# prepare significance computation
#==============================================================================
#%%
res_dic = {"peptide_id":[],
           "exp_values":[],
           "ksuccess":[],
           "ntrials":[],
           "direction":[],
           "pvalue":[],
           "meanlog2fc":[],
           "Sequence": [],
           "ProteinGroup":[],
           "Modifications":[],
           "MH+ [Da]":[],
           "phosphoRS Isoform Probability":[],
           "phosphoRS Site Probabilities":[],
           "seq_window":[],
           "inPSP_sites":[]}

ii = 0
for id_i in all_ids:

    ii += 1
    exp_values_tmp = []
    #init values
    ns1 = np.array([0, 0, 0])
    ns2 = np.array([0, 0, 0])
    ntrials = 0
    if id_i in df_phos1.index:
        # if there is only one value don't do anything
        exp_value1 = avg_entries(df_phos1.loc[id_i])
        #count trials
        ntrials += 1
        ns1 = binom_count(exp_value1, thresholds[0]["upper"])
        exp_values_tmp.append(exp_value1)
        if len(df_phos1.loc[id_i].shape) == 1:
            sequence = df_phos1.loc[id_i]["Sequence"]
            group = df_phos1.loc[id_i]["Protein Group Accessions"]
            mods = df_phos1.loc[id_i]["Modifications"]
            mh = df_phos1.loc[id_i]["MH+ [Da]"]
            isoform = df_phos1.loc[id_i]["phosphoRS Isoform Probability"]
            prob = df_phos1.loc[id_i]["phosphoRS Site Probabilities"]
            window = df_phos1.loc[id_i]["seq_window"]
            inpsps = df_phos1.loc[id_i]["inPSP_sites"]
        else:
            sequence = df_phos1.loc[id_i]["Sequence"][0]
            group = df_phos1.loc[id_i]["Protein Group Accessions"][0]
            mods = df_phos1.loc[id_i]["Modifications"][0]
            mh = df_phos1.loc[id_i]["MH+ [Da]"][0]
            isoform = df_phos1.loc[id_i]["phosphoRS Isoform Probability"][0]
            prob = df_phos1.loc[id_i]["phosphoRS Site Probabilities"][0]
            window = df_phos1.loc[id_i]["seq_window"][0]
            inpsps = df_phos1.loc[id_i]["inPSP_sites"][0]
    else:
        pass

    #check other replicate
    if id_i in df_phos2.index:
        # if there is only one value don't do anything
        exp_value2 = avg_entries(df_phos2.loc[id_i])
        ns2 = binom_count(exp_value2, thresholds[1]["upper"])
        #count trials
        ntrials += 1
        exp_values_tmp.append(exp_value2)
        if len(df_phos2.loc[id_i].shape) == 1:
            sequence = df_phos2.loc[id_i]["Sequence"]
            group = df_phos2.loc[id_i]["Protein Group Accessions"]
            mods = df_phos2.loc[id_i]["Modifications"]
            mh = df_phos2.loc[id_i]["MH+ [Da]"]
            isoform = df_phos2.loc[id_i]["phosphoRS Isoform Probability"]
            prob = df_phos2.loc[id_i]["phosphoRS Site Probabilities"]
            window = df_phos2.loc[id_i]["seq_window"]
            inpsps = df_phos2.loc[id_i]["inPSP_sites"]
        else:
            sequence = df_phos2.loc[id_i]["Sequence"][0]
            group = df_phos2.loc[id_i]["Protein Group Accessions"][0]
            mods = df_phos2.loc[id_i]["Modifications"][0]
            mh = df_phos2.loc[id_i]["MH+ [Da]"][0]
            isoform = df_phos2.loc[id_i]["phosphoRS Isoform Probability"][0]
            prob = df_phos2.loc[id_i]["phosphoRS Site Probabilities"][0]
            window = df_phos2.loc[id_i]["seq_window"][0]
            inpsps = df_phos2.loc[id_i]["inPSP_sites"][0]
    else:
        pass


    ns_12 = ns1 + ns2
    direction = binom_direction(ns_12)

    #compute pvalue
    #k success
    if direction == "neutral":
        ksuccess = 0
    else:
        ksuccess = ns_12.loc[direction]

    x = np.arange(ksuccess, ntrials + 1)
    pvalue = np.sum(stat.binom.pmf(x, ntrials, alpha))
#    print exp_value1, thresholds[0]["upper"]
#    print exp_value2, thresholds[1]["upper"]
#    print ksuccess, "/", ntrials
#    print direction
#    print pvalue
    res_dic["exp_values"].append(exp_values_tmp)
    res_dic["ksuccess"].append(ksuccess)
    res_dic["ntrials"].append(ntrials)
    res_dic["direction"].append(direction)
    res_dic["pvalue"].append(pvalue)
    res_dic["peptide_id"].append(id_i)
    res_dic["meanlog2fc"].append(np.mean(exp_values_tmp))
    res_dic["Sequence"].append(sequence)
    res_dic["ProteinGroup"].append(group)
    res_dic["Modifications"].append(mods)
    res_dic["MH+ [Da]"].append(mh)
    res_dic["phosphoRS Isoform Probability"].append(isoform)
    res_dic["phosphoRS Site Probabilities"].append(prob)
    res_dic["seq_window"].append(window)
    res_dic["inPSP_sites"].append(inpsps)

    if ii == -1:
        break





res_df = pd.DataFrame(res_dic)
res_df["bhd_pvalue"] = multipletests(res_df["pvalue"], method="fdr_bh")[1]
res_df["sig.binom"] = [True if i <= 0.05 else False for i in res_df["bhd_pvalue"]]
res_df = res_df[[u'ProteinGroup', u'Sequence', u'seq_window', u'ksuccess',
                 u'ntrials', u'pvalue', u'bhd_pvalue', u'direction', u'sig.binom',
                 u'meanlog2fc', u'inPSP_sites', u'exp_values', u'peptide_id',
                 u'phosphoRS Isoform Probability',
                 u'phosphoRS Site Probabilities',
                  u'MH+ [Da]', u'Modifications']]
#multiple testing

up_sig = res_df[(res_df["bhd_pvalue"] <= 0.05) & (res_df["direction"] == "up")]
down_sig = res_df[(res_df["bhd_pvalue"] <= 0.05) & (res_df["direction"] == "down")]
neutral_sig = res_df[(res_df["bhd_pvalue"] >= 0.05)]

print "Up", up_sig.shape[0]
print "Down", down_sig.shape[0]
print "Neutral", neutral_sig.shape[0]
#%%
print res_df.head()
res_df.to_csv("/home/sven/data/nurhan/oezge/MayNewData.csv", sep="\t")
test_df = res_df.head.copy()

#%%
y = -np.log10(res_df["bhd_pvalue"] + np.random.normal(0, 0.01, res_df["bhd_pvalue"].shape[0]))
x = res_df["meanlog2fc"]
plt.scatter(x, y, facecolors="none")
plt.axhline(-np.log10(0.05), color="black", lw=1.5)
plt.xlabel("log2fc")
plt.ylabel("-10pvalue")
#%%
np.bincount(res_df["ksucces"])




#==============================================================================
# new stuff
#==============================================================================
def binom_count(exp_value, thresh):
    """
    Return direction and counts for binomial...

    Example:
    -------------------
    exp_value = -1.9
    thresh = 1.01
    """
    ns = pd.Series([0, 0, 0])
    ns.index = ["down", "up", "neutral"]


    if exp_value >= thresh:
        ns["up"] = 1

    elif exp_value <= -thresh:
        ns["down"] = 1

    else:
        ns["neutral"] = 1
    return(ns)


t1 = pd.Series([1, 1, 0])
t1.index = ["down", "up", "neutral"]
t2 = pd.Series([0, 1, 0])
t2.index = ["down", "up", "neutral"]
t3 = pd.Series([1, 0, 1])
t3.index = ["down", "up", "neutral"]
t4 = pd.Series([0, 1, 1])
t4.index = ["down", "up", "neutral"]
binom_direction(t1)
binom_direction(t2)
binom_direction(t3)
binom_direction(t4)

def binom_direction(ns):
    """
    Gets a count array and decides which direction the proteine xpression
    goes.

    Parameters:
    ----------------------
    ns: arr,
         arr of floats...
    """
    if ns.loc["up"] >= ns.loc["down"]:
        return("up")

    elif ns.loc["up"] < ns.loc["down"]:
        return("down")

    else:
        return("neutral")

#    max_idx = np.argmax(ns)
#
#    if max_idx == 0:
#        direction = "up"
#
#    elif max_idx == 1:
#        direction = "down"
#
#    else:
#        direction = "neutral"
#    return(direction)



def avg_entries(serOrdf):
    """
    Gives back a mean or a single value depending on the number of entries
    a specific phosphopeptide has. This is necessary since some peptide seqs
    with the very same modification occur mnultiple times in the output

    Parameter:
    -------------------------------
    serOrdf: series or dataframe,
              series or dataframe from PD. Only the norm_log2HL column is used
    """
    if len(serOrdf) == 1:
        expvalue = serOrdf["norm_log2HL"]

    else:
        #for two values cmpute the mean
        expvalue = serOrdf["norm_log2HL"].mean()

    return(expvalue)


def get_phospho(modsting):
    """
    Function to extract the phosphorylation position from the modification
    string provided by proteomediscoverer. This is necessary to distinguish
    the same peptide sequences but with different mods.

    Parameters:
    ---------------------------------------------
    modstring: str,
               modifcation string from PD

    Example:
    ----------------------------------------------
    tstr = 'N-Term(Dimethyl); K1(Dimethyl); T31(Phospho); K8(Dimethyl); S31(Phospho);'
    get_phospho(tstr)
    >> T31_S31
    """
    p = "([STY]\d+)\(Phospho\)"
    ident = "_".join(re.findall(p, modsting))
    return(ident)


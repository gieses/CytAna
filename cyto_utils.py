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


def venn_diagram(sets, names):
    """
    Plot a Venndiagram

    Parameters:
    -----------------------------
    sets: list of sets,
          elements that should be compared between samples

    names: tuple of str for the sets,
           elements that should be compared between samples
    """
    if len(sets) == 2:
        from matplotlib_venn import venn2, venn2_circles
        f = venn2(sets, names)
        f = venn2_circles(sets)
    elif len(sets) == 3:
        from matplotlib_venn import venn3, venn3_circles
        f = venn3(sets, names)
        f = venn3_circles(sets)
    return (f)

def save_fig(f, out):
    """
    Codes for saving a figre
    """
    f.savefig(out+".png", bbox_inches='tight', pad_inches=0.1)
    f.savefig(out+".svg", bbox_inches='tight', pad_inches=0.1)
    f.savefig(out+".pdf", bbox_inches='tight', pad_inches=0.1)
    f.clf()
    plt.close()

def get_color_map(name, map_type, number, reverse=False, get_map=False):
    """
    Parameters
    ----------
    name : str
        Name of color map. Use `print_maps` to see available color maps.
    map_type : {'Sequential', 'Diverging', 'Qualitative'}
        Select color map type.
    number : int
        Number of defined colors in color map.
    reverse : bool, optional
        Set to True to get the reversed color map.

    All colors:
    Sequential
    Blues     :  {3, 4, 5, 6, 7, 8, 9}
    BuGn      :  {3, 4, 5, 6, 7, 8, 9}
    BuPu      :  {3, 4, 5, 6, 7, 8, 9}
    GnBu      :  {3, 4, 5, 6, 7, 8, 9}
    Greens    :  {3, 4, 5, 6, 7, 8, 9}
    Greys     :  {3, 4, 5, 6, 7, 8, 9}
    OrRd      :  {3, 4, 5, 6, 7, 8, 9}
    Oranges   :  {3, 4, 5, 6, 7, 8, 9}
    PuBu      :  {3, 4, 5, 6, 7, 8, 9}
    PuBuGn    :  {3, 4, 5, 6, 7, 8, 9}
    PuRd      :  {3, 4, 5, 6, 7, 8, 9}
    Purples   :  {3, 4, 5, 6, 7, 8, 9}
    RdPu      :  {3, 4, 5, 6, 7, 8, 9}
    Reds      :  {3, 4, 5, 6, 7, 8, 9}
    YlGn      :  {3, 4, 5, 6, 7, 8, 9}
    YlGnBu    :  {3, 4, 5, 6, 7, 8, 9}
    YlOrBr    :  {3, 4, 5, 6, 7, 8, 9}
    YlOrRd    :  {3, 4, 5, 6, 7, 8, 9}
    Diverging
    BrBG      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    PRGn      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    PiYG      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    PuOr      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    RdBu      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    RdGy      :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    RdYlBu    :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    RdYlGn    :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    Spectral  :  {3, 4, 5, 6, 7, 8, 9, 10, 11}
    Qualitative
    Accent    :  {3, 4, 5, 6, 7, 8}
    Dark2     :  {3, 4, 5, 6, 7, 8}
    Paired    :  {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}
    Pastel1   :  {3, 4, 5, 6, 7, 8, 9}
    Pastel2   :  {3, 4, 5, 6, 7, 8}
    Set1      :  {3, 4, 5, 6, 7, 8, 9}
    Set2      :  {3, 4, 5, 6, 7, 8}
    Set3      :  {3, 4, 5, 6, 7, 8, 9, 10, 11, 12}

    """
    if get_map is True:
        return(brewer2mpl.get_map(name, map_type, number).mpl_colormap)
    else:
        return(brewer2mpl.get_map(name, map_type, number).mpl_colors)


def read_kinase(file_loc="E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\Phosphosite_data\\Kinase_Substrate_Dataset"):
    """
    Readsa data from PhosphositePLus and offers it in a dictionary format
    """
    df = pd.read_csv(file_loc, sep="\t")
    return (df)

def read_phosphosite(file_loc="E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\Phosphosite_data\\Kinase_Substrate_Dataset"):
    """
    Readsa data from PhosphositePLus and offers it in a dictionary format
    """
    df = pd.read_csv(file_loc, sep="\t")
    df = df[df["ORGANISM"] == "human"]
    df = df.set_index("ACC_ID")
    return (df)


def read_kinome_explorer(infile):
    """
    Reads the kinome explorer file
    """

    df_kinome = pd.read_csv(infile, sep="\t")
    df_kinome = df_kinome.set_index("substrate")
    return(df_kinome)

def add_known_pssite(phospho_df, phosphositeDB, fasta_dic):
    """
    Checks if the presented phosphopeptide is already known in the
    PSP database.

    Parameters:
    ----------------------
    phospho_df: str,
             uniprot accession
    phosphositeDB: dataframe,
           indexed dataframe containing the PSP data
    """
    index_dic = {i for i in phosphositeDB.index}
    uniprot_kinome = []
    position_kinome = []
    aminoacid_kinome = []
    global_bools = []
    all_bools = []
    all_phospho = []
    for row in phospho_df.iterrows():
        temp_bool = []
        temp_phospho = []
        #get current phosphorylation data
        accession = row[1]["Protein Group Accessions"]
        #check if accession is phosphosite database
        if (accession in index_dic) and (accession in fasta_dic):
            phospho_peptide =  row[1]["Sequence"]
            # get modification position (in peptide sequence) and amino acid
            sites =  get_phospho(row[1]["Modifications"], position=True)
            sites_aminoacid = get_phospho(row[1]["Modifications"], position=False)

            #for eac phosphorylation site get the KinomeXplorer id string
            for phossite, phosaa in zip(sites.split(";"), sites_aminoacid.split(";")):

                kinome_t = fasta_dic[accession].replace("I", "L").index(
                         phospho_peptide.replace("I", "L")) + int(phossite)# - 1
                #sequence = fasta_dic[accession]
                #print accession, phospho_peptide, phossite,  sequence[kinome_t]
                uniprot_kinome.append(accession)
                position_kinome.append(kinome_t)
                aminoacid_kinome.append(phosaa)
                #create an identifier and check if the phossite
                #is present in the database
                identifier_PSPS = "{}{}-p".format(phosaa, kinome_t)
                temp_phospho.append(identifier_PSPS)
                #needed to circument pandas datafram behavior (if single
                #entry only a string is returned)
                try:
                    identifier_temp = phosphositeDB.ix[accession].values
                except:
                    identifier_temp = [phosphositeDB.ix[accession]]

                if identifier_PSPS in identifier_temp:
                    temp_bool.append(True)

                else:
                    temp_bool.append(False)


            #global bools indicate if there was support for any of the sites
            #in PSPS. True can come from a single "verified" site
            if np.any(temp_bool):
                global_bools.append(True)
            else:
                global_bools.append(False)
            all_phospho.append(";".join(temp_phospho))
            all_bools.append(";".join([str(booli) for booli in temp_bool]))

        #if accessions is not known, do nothing
        else:
            all_phospho.append("")
            global_bools.append(False)
            all_bools.append(False)
            continue

    phospho_df["inPSP_any"] = global_bools
    phospho_df["inPSP_individual"] = all_bools
    phospho_df["inPSP_sites"] = all_phospho

def mad(arr, correction=1.482):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation

    Computes the MAD of a given array
    Paramters:
    --------------
    arr: array-like,
         numbers

    correction: float,
                correction factor for estimating normal parameters

    Note:
    ------------------
    When normal distribution paramters are estimated with the mad, a correction
    factor (1.482) needs to be applied to match the estimated distribution
    mean
    """
    # should be faster to not use masked arrays.
    arr = np.ma.array(arr).compressed()
    med = np.median(arr)
    return np.median(np.abs(arr - med) * 1.482)


def estimate_normal_params(x, outlier=True, value=3.):
    """
    Estimate the mean (as median) and veriance (as MAD)

    Paramters:
    ----------------------

    x: array,
       float
     outlier: bool,
              with or without outlier removal
    value: float,
           factor for outlier removal (value * sd +. median is removed)

    Returns:
    ---------------------
    (loca, scale): float, float
                  tuple of loc and scale estimation
    """
    if outlier:
        mu, sd = norm.fit(x)
        mu = np.median(x)
        sd = mad(x)
        t = x
        t = t[(t >= mu - value * sd) & (t <= mu + value * sd)]
        loc = np.mean(t)
        scale = np.std(t)
    else:
        loc = np.median(x)
        scale = mad(x)

    return (loc, scale)

def get_fasta_dic(path="/home/sven/Dropbox/VX data for Sven/sprot_2014_08_2014_11.fasta"):
    """
    Retrieve a <accesion>:<sequence> mapping from a fasta file.
    """
    try:
        seq_dic = {}
        fasta = HTSeq.FastaReader(path)
        for fentry in fasta:
            seq_dic[fentry.name.split("|")[1]] = fentry.seq
    except:
        seq_dic = {}
        fasta = HTSeq.FastaReader("E:\\cloud_space\\Dropbox\\VX data for Sven\\sprot_2014_08_2014_11.fasta")
        for fentry in fasta:
            seq_dic[fentry.name.split("|")[1]] = fentry.seq
    return(seq_dic)


def get_kinomexplorer_output(phospho_df, fasta_dic, outfile, fasta_file):
    """
    Creates an outputfie for the kinome explorer

    Parameter:
    -----------------------
    phospho_df: df,
                 PD dataframe
    fasta_dic: dic,
               <uniprot>:<sequence> dictionary
    outfile: str,
             destination output

    Output:
    -----------------------
    writes an outputfile for KinaseXplorer
    Accession    position    aminoacid(phospho)

    """
    stats_missed = 0
    c = 0
    #iterate over dataframe


    uniprot_kinome = []
    position_kinome = []
    aminoacid_kinome = []
    for row in phospho_df.iterrows():
        #get current phosphorylation data
        accession = row[1]["Protein Group Accessions"]
        phospho_peptide =  row[1]["Sequence"]
        # get modification position (in peptide sequence) and amino acid
        sites =  get_phospho(row[1]["Modifications"], position=True)
        sites_aminoacid = get_phospho(row[1]["Modifications"], position=False)
        if accession in fasta_dic:
            #for eac phosphorylation site get the KinomeXplorer id string
            for phossite, phosaa in zip(sites.split(";"), sites_aminoacid.split(";")):

                kinome_t = fasta_dic[accession].replace("I", "L").index(
                         phospho_peptide.replace("I", "L")) + int(phossite)# - 1
                #sequence = fasta_dic[accession]
                #print accession, phospho_peptide, phossite,  sequence[kinome_t]
                uniprot_kinome.append(accession)
                position_kinome.append(kinome_t)
                aminoacid_kinome.append(phosaa)
        else:
            stats_missed +=1
        c += 1

    #transform to handy dataframe
    kinome_df = pd.DataFrame()
    kinome_df["Accesion"] = uniprot_kinome
    kinome_df["Position"] = position_kinome
    kinome_df["PhosAA"] = aminoacid_kinome
    print "{} sequences were not found in the database".format(stats_missed)
    #kinome_df.to_clipboard(index=False)
    print "write to file: {}".format(outfile)
    kinome_df.to_csv(outfile, sep="\t", index=False)

    seq_dic = {}
    fasta = HTSeq.FastaReader(fasta_file)
    for fentry in fasta:
        seq_dic[fentry.name.split("|")[1]] = fentry.seq
    return(seq_dic)



def get_sequence_window(peptide, accession, fastadb, sites, window=7,
                        verbose=False):
    """
    Function to retrieve a sequence window with +- window aminoa cids from the
    phosphorylated amino acid


    Paramters:
    -------------------------------------
    peptide: str,
             peptide sequence
    accesion: str,
              proteina ccesion
    fastadb: dic,
             <accession>:sequence dictionary
    window: int,
            number of amino acids before / after phosphorylated residue
    """

    #    peptide = "EVYELLDSPGK"
    #    accession = "P22234"
    #    phossite = 8
    #
    #    peptide ="MATAEVLNIGKK"
    #    phossite = 1 #5
    #
    #    peptide ="QADKKIRECNL"
    #    phossite = 11 #
    peptide = peptide.replace("I", "L")
    window_seqs = []
    for site in sites.split(";"):
        phossite = int(site)
        if accession in fastadb:
            center = fastadb[accession].replace("I", "L").index(peptide) + phossite
            #need some clipping function because the sequences are "__" if
            #protein
            # is 'shorter'
            if center < window:
                prefix = "_" * (window - center + 1) + fastadb[accession][:center]
                suffix = fastadb[accession][center: center + window]
                window_seqs.append(prefix + suffix)

            elif center > len(fastadb[accession]) - window:
                diff = np.abs((len(fastadb[accession]) - window) - center)
                prefix = fastadb[accession][center - window - 1: center]
                suffix = fastadb[accession][center:] + "_" * (diff)
                window_seqs.append(prefix + suffix)
            else:
                # all fine, just get the window
                window_seqs.append(fastadb[accession][center - 1 - window: center + window])
            if verbose:
                    print (window_seqs)
                    print (len(window_seqs))
                    print fastadb[accession]
                    print fastadb[accession][center]
        else:
            window_seqs = ["NA"]

    return(";".join(window_seqs))


def explorative_scatter(phospho_peptides):
    f, ax = plt.subplots(1, 2, figsize=(18, 12))
    ax[0].scatter(phospho_peptides["log10HL"],
                  phospho_peptides["log10ML"], c="g", alpha=0.3,
                  label="MitCyt/IntCyt")

    ax[0].scatter(phospho_peptides["log10HL"],
                  phospho_peptides["log10HM"], c="b", alpha=0.3,
                  label="MitCyt/MitInt")
    ax[0].legend(loc="lower left")
    ax[0].set(xlabel="log10 ratio (HL)", ylabel="log10 ratio (ML and HM)")

    ax[1].scatter(phospho_peptides["log10ML"],
                  phospho_peptides["log10HM"], c="y", alpha=0.3,
                  label="IntCyt/MitInt")

    ax[1].scatter(phospho_peptides["log10ML"],
                  phospho_peptides["log10HL"],
                  c="g", alpha=0.3, label="IntCyt/MitCyt")

    ax[1].legend(loc="lower left")
    ax[1].set(xlabel="log10 ratio (ML)", ylabel="log10 ratio (HM and HL)")
    #plt.savefig("E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\exploratory_figures_tables\\{}_compare_states.png".format(filename))


def explorative_violin(t):
    bmap = brewer2mpl.get_map('Paired', 'Qualitative',8).mpl_colors
    plt.figure(figsize=(18, 12))
    sns.boxplot(data=t)
    plt.axhline(0, alpha=0.6, c="k")
    plt.axhline(0.5, alpha=0.6, c="k")
    plt.xticks(np.arange(1, len(t)+1), xlabels)
    plt.xlabel("experiment type (peptide specific)")
    plt.ylabel("log10(foldchange)")
    plt.title("unnormalized ratios")
    #plt.savefig("E:\\cloud_space\\Dropbox\\VX data for Sven\\SVEN\\exploratory_figures_tables\\{}_visualize_expression_violin.png".format(filename))


def get_non_null(dataframe, column=["log10HL", "log10ML", "log10HM"]):
    """
    Gets all the values from a pandas series that are not null.

    Parameters:
    -------------------------------
    dataframe: pd.df
               dataframe of peptide / phosphopeptides
    """
    return(dataframe[pd.notnull(dataframe[column].values)][column].values)


def normalize_df(phospho_df, reference, columns, verbose=False):
    """
    Normalizes a dataframe expression column

    Parameters:
    -----------------------
    phospho_df: df
                phosphodata dataframe
    referene: df,
              non-phosphorylated reference
    columns: list,
             string of columns
    """

    for ratio in columns:
        reference_med = np.median(reference[ratio])
        print "shift distribution by: {}".format(reference_med)
        phospho_df["norm_"+ratio] = phospho_df[ratio] - reference_med


def get_phospho(mod_str, position=False):
    """
    Extracts the phospho identifier from a modification string from PD input

    Parameter:
    -------------------
    mod_str: str,
             input parameter from PD output ("ModificatioN" column)

    position: bool,
              1 - get the position
              0 - get the amino acid
    """
    # example
    #modstr="N-Term(Dimethyl:2H(6)13C(2)); T11(Phospho); K13(Dimethyl:2H(6)13C(2))"
    global multiple_sites
    if position:
        phos_type = re.findall("[STY](\d+)\(Phospho\)", mod_str)
        return(";".join([i for i in phos_type]))
    else:
        phos_type = re.findall("([STY])\d+\(Phospho\)", mod_str)
    if len(phos_type) > 1:
        multiple_sites += 1
        phos_type = ";".join([i[0] for i in phos_type])
    else:
        phos_type = phos_type[0]

    return (phos_type)


def filters():
    print "Get filtered peptides...."
    #False if Nan
    HL_filter = pd.notnull(peptides_df["Heavy/Light"])
    ML_filter = pd.notnull(peptides_df["Medium/Light"])

    filter_ = [True if i == True or j == True else False for i, j in zip(HL_filter, ML_filter)]
    filtered_peptides_df = peptides_df[filter_].copy()
    filtered_peptides_df["Heavy/Medium"] = filtered_peptides_df["Heavy/Light"] / filtered_peptides_df["Medium/Light"]
    filtered_peptides_df["isPhospho"] = pd.notnull(filtered_peptides_df["phosphoRS Isoform Probability"])


    log10fphos = [pd.notnull(i) and pd.notnull(j) and pd.notnull(k) for i, j, k
                   in
                   zip(phospho_peptides["log10HL"].values,
                       phospho_peptides["log10ML"].values,
                       phospho_peptides["log10HM"].values)]

    log10freg = [pd.notnull(i) and pd.notnull(j) and pd.notnull(k) for i, j, k
                   in
                   zip(regular_peptides["log10HL"].values,
                       regular_peptides["log10ML"].values,
                       regular_peptides["log10HM"].values)]

    filtered_regular = regular_peptides[log10freg]
    filtered_phosphor = phospho_peptides[log10fphos]

    imputer = Imputer(missing_values='NaN', strategy="knn",
                      axis=0, n_neighbors=5)
    X_impute = imputer.fit(peptides_df[ratios_default]).transform(peptides_df[ratios_default])


def computelog(dataframe, ratio_columns=["Heavy/Light", "Medium/Light",
                                         "Heavy/Medium"]):
    """
    Comute the log10 columns in the dataframe.

    Parameters:
    -----------------------------------------------
    dataframe: df,
               pandas dataframe with ratio columns
    ratio_columns: list,
                    list of ratio names, format important. For example:
                    ["Heavy/Light", "Medium/Light", "Heavy/Medium"] works

    """
    for rc in ratio_columns:
        short = rc[0] + rc.split("/")[1][0]
        dataframe["log10{}".format(short)] = np.log10(dataframe[rc])


def get_nonNan(x):
    """
    Return nonNan values
    """
    return(x[pd.notnull(x)])


def plt_hist(xin, outfile, column):
    """
    Plot the histogram
    """
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    loc, scale = estimate_normal_params(xin, outlier=False)
    ax.hist(xin, bins=40, normed=True, alpha=0.7)
    ax.set_xlim(-3, 3)
    xmin, xmax = ax.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, loc, scale)
    ax.plot(x, p, 'k-', linewidth=2, alpha=0.8)
    ax.set_xlabel("log2 (foldchange)")
    ax.set_ylabel("Density")
    sns.despine()
    ax.axvline(loc + scale, ls="--", lw=2, color="k", alpha=0.7)
    ax.axvline(loc - scale, ls="--", lw=2, color="k", alpha=0.7)
    ax.axvline(loc, ls="--", lw=2, color="k", alpha=0.7)
    save_fig(f, ("{}_{}".format(outfile, column)))


def plt_correlation(x, y, outfile):
    xval = (x + y) / 2
    yval = x - y
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    ax.plot(xval, yval, 'ko', linewidth=2, alpha=0.8)
    ax.set_xlabel(" avg(Light + Heavy)")
    ax.set_ylabel("heavy - light")
    sns.despine()
    save_fig(f, ("{}_{}".format(outfile)))


def analyze_ratio(phospho_peptides, regular_peptides, column, alpha,
                  outfile, onlysig=False):
    """
    Routine for analyzing a specific ratio.

    Parameters:
    -------------------------------
    phospho_peps: df,
                  phosphopeptide dataframe

    regular_peps: df,
                  regular peptides dataframe

    column: str,
                  column that is used for the analysis

    alpha: float,
           error level

    onlysig: bool,
             if False report only significant peptides to outputfile
    """
    phospho_filtered = phospho_peptides[pd.notnull(phospho_peptides[column])].copy()
    regular_filtered = regular_peptides[pd.notnull(regular_peptides[column])].copy()

    plt_hist(phospho_filtered[column].values, outfile+"_phospho", column)
    plt_hist(regular_filtered[column].values, outfile+"_regular", column)


    normalize_df(phospho_filtered, regular_filtered, [column])
    normalize_df(regular_filtered, regular_filtered, [column])

    loc, scale = estimate_normal_params(regular_filtered["norm_"+column], False)
    lower_bound = stat.norm.ppf(alpha/2, loc=loc, scale=scale)
    upper_bound = stat.norm.ppf(1 - alpha/2, loc=loc, scale=scale)
 #%%
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    temp_df = pd.DataFrame()
    temp_df["Log10Ratio"] = [value_i for value_i in list(phospho_filtered[column])
                                            + list(phospho_filtered["norm_"+column])
                                            + list(regular_filtered[column])
                                            + list(regular_filtered["norm_"+column])]

    temp_df["identifier"] = \
    ["phospho_filtered" for value_i in phospho_filtered[column]] \
    + ["phospho_filtered_norm" for value_i in phospho_filtered["norm_"+column]] \
    + ["regular_filtered" for value_i in regular_filtered[column]] \
    + ["regular_filtered_norm" for value_i in regular_filtered["norm_"+column]]


    #flierprops = dict(marker='o', markersize=5, color="k")
    ax = sns.boxplot(x="Log10Ratio", y="identifier", data=temp_df)

#    ax.boxplot([phospho_filtered[column], phospho_filtered["norm_"+column],
#                 regular_filtered[column], regular_filtered["norm_"+column]])
    ax.axvline(0, color="k", ls="--", alpha=0.8)
    ax.set_xlabel("norm "+ column)
    ax.set_ylabel("peptide type")
    ax.axvline(lower_bound, color="k", ls="--", alpha=0.8)
    ax.axvline(upper_bound, color="k", ls="--", alpha=0.8)
    ax.set_yticks([0, 1, 2, 3],)
    ax.set_yticklabels( ["Phospho", "Norm Phospho", "Regular", "Norm Regular"])
    sns.despine()

    save_fig(f, ("{}_{}_boxplot".format(outfile, column)))
     #%%

    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    ax.scatter(phospho_filtered["norm_"+column], phospho_filtered["PEP"])
    ax.set_xlabel("norm_"+column)
    ax.set_ylabel("PEP score")
    #phospho_filtered["Heavy/Light Variability [%]"])
    ax.axvline(lower_bound, c="k", ls="--")
    ax.axvline(upper_bound, c="k", ls="--")
    ax.axhline(alpha, c="k", ls="--")
    sns.despine()
     #%%
    save_fig(f, "{}_{}_scatter".format(outfile, column))


    #get annotation for significant peptides
    phospho_filtered["significant"] = [True if (i >= upper_bound or i <= lower_bound)
                                       else False for i in
                                       phospho_filtered["norm_"+column]]

    #get annotation for the direction, i.e. is the peptide a candidate
    # for up (fc >= 0) or down (fc <0) regulation
    phospho_filtered["direction"] = ["Up" if i >= 0 else "Down" for
                                     i in phospho_filtered["norm_"+column]]

    if onlysig:
        significants = phospho_filtered[phospho_filtered["significant"] == True]
    else:
        significants = phospho_filtered.copy()
        significants.sort_values(by="significant", inplace=True)
 #%%
    #significants.to_csv("{}_{}.csv".format(outfile, column), sep="\t")
    return(significants)

def add_predicted_groups(phospho_peptides, kinome_df):
    """
    Adds motif predictions from the kinase explorer
    """
    groups = []
    idx_keys = {i for i in kinome_df.index}
    for row_i in phospho_peptides.iterrows():
        accession = row_i[1]["Protein Group Accessions"]
        pos = row_i[1]["inPSP_sites"].split("-")[0][1:]
        if accession in idx_keys:
            tgroups = kinome_df.ix[accession][kinome_df.ix[accession]["position"]]
        else:
            pass

def add_seq_window(phospho_peptides, FASTA_dic):
    """
    Function that adds a sequence window column
    """
    phospho_peptides["seq_window"] = [get_sequence_window(pep, acc, FASTA_dic,
                                 get_phospho(site, position=True))
                                 for pep, acc, site in zip(
                                 phospho_peptides["Sequence"],
                                 phospho_peptides["Protein Group Accessions"],
                                 phospho_peptides["Modifications"])]


def get_phosphorylated_stats(phospho_peptides):
    """
    Get statistics about phosphopeptide occurrences

    Parameters:
    ------------------------
    phospho_peptides: df
                      PD dataframe
    """
    print "Get Phosphostats peptides...."
    phospho_peptides["PhosAA"] = [get_phospho(modi) for modi in
                                  phospho_peptides["Modifications"]]
    print "Multiple sites: {}".format(multiple_sites)
    frequency_tab = np.round(phospho_peptides["PhosAA"].value_counts(
                             normalize=True), 3)
    return (frequency_tab)


def probe_protein_abundance(phos_df, reg_df, prot_df, column):
    """
    **experimental only - not used**
    """
    intersection = np.intersect1d(phospho_peptides["Sequence"],
                                  regular_peptides["Sequence"])
    intersection_val = float(len(intersection)) / phospho_peptides.shape[0]
    print "{}% (#{})".format(np.round(intersection_val, 2), len(intersection))

    new_prot_df = proteins_df.set_index("Accession").copy()
    new_phos_df = phospho_peptides.copy()
    new_reg_df = regular_peptides.copy()
    new_phos_df = new_phos_df.set_index(["Sequence"])
    new_reg_df = new_reg_df.set_index(["Sequence"])

    new_phos_df = new_phos_df.loc[intersection]

    phosphopep = []
    regulpep = []
    protein_ratio = []
    for row in new_phos_df.iterrows():
        try:
            accessions = row[1]["Protein Group Accessions"].split(";")
            prot = new_prot_df.loc[accessions[0]]
            regular_p = new_reg_df.loc[row[0]]
            t = [row[1][column], prot[column], regular_p[column]]
            if sum(pd.isnull(t)) == 0:
                phosphopep.append(row[1][column])
                regulpep.append(regular_p[column])
                protein_ratio.append(prot[column])
            else:
                continue
        except:
            pass

    x = np.array(phosphopep)
    y = np.array(regulpep)
    z = np.array(protein_ratio)

    Lp = (z - y) / (x - z)
    Hp = (x * (z - y)) / (y * (x - z))
    L = Lp / (1 + Lp)
    H = Hp / (1 + Hp)

    plt.plot(phosphopep, regulpep, 'ro')
    plt.xlabel("phospho")
    plt.ylabel("regular")
    plt.plot(phosphopep, protein_ratio, 'bo')
    plt.xlabel("phospho")
    plt.ylabel("protein")
    plt.title("Correlation: Phosphopeptide vs. Protein ratio")

    plt.plot(regulpep, protein_ratio, 'go')
    plt.xlabel("regular")
    plt.ylabel("protein")
    from scipy.stats import pearsonr
    print "Phospho vs. Regular", pearsonr(phosphopep, regulpep)
    print "Phospho vs. Protein", pearsonr(phosphopep, protein_ratio)
    print "Regular vs. Protein", pearsonr(regulpep, protein_ratio)
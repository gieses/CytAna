# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 22:33:32 2016

@author: Hanjo
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import cyto_utils as cutils
from matplotlib.ticker import MaxNLocator
from scipy.stats import norm
from scipy import stats
import brewer2mpl
import matplotlib as mpl
fontsize = 26
JR_style = {'axes.axisbelow': True,
             'axes.grid': False,
             'axes.labelcolor': '.15',
             'axes.linewidth': 1,
             'font.family': ['sans-serif'],
             'font.sans-serif': ['Arial',
              'Liberation Sans',
              'Bitstream Vera Sans',
              'sans-serif'],
             'image.cmap': 'Greys',
             'legend.frameon': False,
             'legend.numpoints': 1,
             'legend.scatterpoints': 1,
             'lines.solid_capstyle': 'round',
             'xtick.direction': 'out',
             'xtick.major.size': 8,
             'xtick.minor.size': 5,
             'ytick.direction': 'out',
             'ytick.major.size': 8,
             'ytick.minor.size': 5,
             'legend.fontsize': 20,
             'font.size': fontsize,
             'axes.titlesize': fontsize-5,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'axes.labelsize': fontsize,
              'svg.fonttype': 'none'}
sns.set_context("notebook")
sns.set_style("ticks")
mpl.rcParams.update(JR_style)
#%%
#sns.set_style("ticks")
#sns.set_context("talk", font_scale=1.8)

def fig_ratiohisto(xvalues, column, figpath):
    #%%
    alpha = 0.05
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    ax.hist(xvalues, bins=80, normed=True, alpha=0.7)
    maxx = np.max([np.max(xvalues), abs(np.min(xvalues))]) + 0.1
    x = np.linspace(-maxx, maxx, 100)
    loc, scale = cutils.estimate_normal_params(xvalues)
    lower_bound = stats.norm.ppf(alpha/2, loc=loc, scale=scale)
    upper_bound = stats.norm.ppf(1 - alpha/2, loc=loc, scale=scale)
    p = norm.pdf(x, loc, scale)
    ax.plot(x, p, '--', linewidth=2, alpha=0.8, lw=2, c="k")
    ax.set_xlim(-1.5, 1.5)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.yaxis.set_major_locator(MaxNLocator(4))
    ax.set(xlabel="log10 (foldchange)", ylabel="Density", title="Normalized peptides ratio {}".format(column))
    ax.axvline(lower_bound, ls="--", lw=2, color="red", alpha=0.7)
    ax.axvline(upper_bound, ls="--", lw=2, color="red", alpha=0.7)
    sns.despine()
    cutils.save_fig(f, figpath+"162_FitHistogram_{}".format(column))
    return(lower_bound, upper_bound)
    #%%

def fig_boxplotcomparison(regular, phospho, sheet, column, bounds, figpath):
    """
    Figure 2c from Aurora paper
    """
    #%%
    xvals = []
    yvals = []
    for i, j in zip([regular[sheet], regular[column], phospho[sheet].values, phospho[column].values], ["peptides", "norm_peptides", "phospho", "norm_phospho"]):
        xvals.extend(i)
        yvals.extend([j]*len(i))

    df = pd.DataFrame([xvals, yvals]).transpose()
    df.columns = [column, "peptide type"]
    #%%
    lower_bound = bounds[0]
    upper_bound = bounds[1]
    bmap = brewer2mpl.get_map('Paired', 'Qualitative', 6).mpl_colors
    f, ax = plt.subplots(1, figsize=(11.69, 8.27))
    sns.boxplot(x="peptide type", y=column, data=df, palette=bmap)
    ax.set(xticks=[0, 1, 2, 3],
           xticklabels=["peptides", "peptides \n(normalized)", "phospho-\npeptides", "phosphopeptides \n(normalized)"])
    ax.set(ylim=(-2, 2), ylabel="log10 (fold change)", title=column)
    ax.axhline(lower_bound, ls="--", lw=2, color="red", alpha=0.7)
    ax.axhline(upper_bound, ls="--", lw=2, color="red", alpha=0.7)
    ax.yaxis.set_major_locator(MaxNLocator(4))
    sns.despine()
    cutils.save_fig(f, figpath+"162_BoxplotNorm_{}".format(column))
    #%%

def fig_phosints(infile1, infile2, figpath):
    """
    Visualize the non-phospgorylated vs. phosphorylated peptide ratios.
    """
    #sheetnames = ["log10HL", "log10HM", "log10ML"]
    sheetnames = ["log10HL"]
    for sheet in sheetnames:
        phospho = pd.read_excel(infile1, sheetname=sheet)
        regular = pd.read_excel(infile2, sheetname=sheet)
        column = "norm_"+ sheet
        #%%

        #======================================================================
        #         Scatterplot: General trend
        #======================================================================
        maxx = np.max([np.max(phospho[column]), abs(np.min(phospho[column]))]) + 0.1
        f, ax = plt.subplots(1, figsize=(11.69, 8.27))
        ax.scatter(regular[column], np.log10(regular["Area"]), s=80, lw=1,
                   alpha=0.7, label="Regular peptides")
        ax.scatter(phospho[column], np.log10(phospho["Area"]), c="green",
                   s=80, lw=1, alpha=0.7, label="phospho peptides")
        ax.set(title="Normalized peptide Ratios: {}".format(sheet),
               xlabel="log10 (fold change)", ylabel="log10 (intensity)",
                xlim=(-maxx, maxx))
        ax.legend(loc="upper left")
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        sns.despine()
        cutils.save_fig(f, figpath+"162_Intensityscatter_{}".format(column))
        #%%
        #======================================================================
        #         Histogram
        #======================================================================
        #%%
        xvalues = regular[column].values
        # ratio histogram
        bounds = fig_ratiohisto(xvalues, column, figpath)

        #boxplot over all ratios
        fig_boxplotcomparison(regular, phospho, sheet, column, bounds, figpath)


##==============================================================================
## 162 data
##==============================================================================
#figpath = "D:\\Sven\\Dropbox\\VX data for Sven\\SVEN\\Figures\\CytoNoc\\"
#infile1 = "D:\Sven\Dropbox\VX data for Sven\SVEN\data\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide_overview.xlsx"
#infile2 = "D:\Sven\Dropbox\VX data for Sven\SVEN\data\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide\ALL162files_Mit_Int_CytoPurA_mix_SCX_TiO2_peptide_nonhospho_overview.xlsx"
#fig_phosints(infile1, infile2, figpath)
#
###==============================================================================
### Figure 1: regular peptides vs. phosphopeptides
###==============================================================================
#path = "D:\\Sven\Dropbox\\VX data for Sven\\SVEN\data\\Mit_Int_CytoNoc_mix_IEF_TiO2_peptide\\"
#figpath = "D:\\Sven\\Dropbox\\VX data for Sven\\SVEN\\Figures\\CytoPurA\\"
#infile1 = path + "Mit_Int_CytoNoc_mix_IEF_TiO2_peptide_overview.xlsx"
#infile2 = path + "Mit_Int_CytoNoc_mix_IEF_TiO2_peptide_nonhospho_overview.xlsx"
#fig_phosints(infile1, infile2, figpath)


#==============================================================================
# new data
#==============================================================================
path = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\May_2016_new cytokinesis files\\"
figpath = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\May_2016_new cytokinesis files\\PH_NL"
infile1 = path + "PH_NL_peptides_overview.xlsx"
infile2 = path + "PH_NL_peptides_nonhospho_overview.xlsx"
fig_phosints(infile1, infile2, figpath)


path = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\May_2016_new cytokinesis files\\"
figpath = "D:\\Sven\\Dropbox\\shared_folders\\nurhan\\VX data for Sven\\SVEN\May_2016_new cytokinesis files\\PL_NH"
infile1 = path + "PL_NH_peptides_nonhospho_overview.xlsx"
infile2 = path + "PL_NH_peptides_overview.xlsx"
fig_phosints(infile1, infile2, figpath)



import pandas as pd

conn = "localhost"
dbname ="testdb"
df = pd.read_csv("C:\\temp\\test.csv", sep="\t")
df.to_sql()
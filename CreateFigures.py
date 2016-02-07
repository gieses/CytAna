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
ax.plot()

excel.columns


plt.hist(excel[excel["direction"]=="Down"]["norm_log10HL"].values, bins=20, label="Down")
plt.hist(excel[excel["direction"]=="Up"]["norm_log10HL"].values, bins=20, label="Up")
plt.legend(loc="upper left")

excel2 = excel[excel["significant"] == True]
plt.hist(excel2[excel2["direction"]=="Down"]["norm_log10HL"].values, bins=20, label="Sig. Down")
plt.hist(excel2[excel2["direction"]=="Up"]["norm_log10HL"].values, bins=20, label="Sig. Up")
plt.legend(loc="upper left")

excel2 = excel[excel["significant"] == False]
plt.hist(excel2[excel2["direction"]=="Down"]["norm_log10HL"].values,  label="Sig. Down")
plt.hist(excel2[excel2["direction"]=="Up"]["norm_log10HL"].values,  label="Sig. Up")
plt.xlim(-2, 2)
plt.legend(loc="upper left")


plt.boxplot([phospho_filtered["norm_log10HM"],
             phospho_filtered["log10HM"],
             regular_filtered["log10HM"],
             regular_filtered["norm_log10HM"]])
plt.xticks([1,2,3,4], ["phospho norm", "phospho", "regular norm", "regular"])

plt.hist(phospho_filtered["norm_log10HM"].values, bins=30)
plt.hist(regular_filtered["norm_log10HM"].values, bins=30)
plt.hist(regular_filtered["log10HM"].values, bins=30)
plt.hist(phospho_filtered["norm_log10HM"].values, bins=30)
plt.xlim(-1.2, 1.2)

sns.set_style("ticks")
sns.set_context("talk", font_scale=1.4)
plt.scatter(regular_filtered["norm_log10HM"],
            np.log10(regular_filtered["Area"]), s=80, lw=1, alpha=0.7,
            label="Regular peptides")
plt.scatter(phospho_filtered["norm_log10HM"],
            np.log10(phospho_filtered["Area"]), c="green", s=80, lw=1,
            alpha=0.7, label="phospho peptides")
plt.title("Normalized peptide Ratios: HM")
plt.xlabel("log2 (fold change)")
plt.ylabel("log10 (intensity)")
plt.xlim(-3, 3)
sns.despine()
plt.show()
# %%
import pandas as pd

df = pd.read_csv("Filter_tr_genotypes/trgt_denovo_report.tsv", sep="\t")
# %%
print(df.head())
#df.describe()
# %%
len(df)
# %%
df_filtered_denovo = df[(df["denovo_coverage"] >=3) & (df["child_ratio"] >= 0.2)]
len(df_filtered_denovo)
df_filtered_denovo.head()
df_filtered_denovo = df_filtered_denovo[~df_filtered_denovo["trid"].str.contains("Un|random", regex=True)]
df_filtered_denovo.head()
len(df_filtered_denovo)
df_filtered_denovo[["child_ratio","child_coverage","per_allele_reads_child","allele_ratio","allele_coverage","denovo_coverage",]].head()
# %%

# annotate with depth and filter on depth
df_filtered_denovo["mom_dp"] = df_filtered_denovo["per_allele_reads_mother"].apply(
        lambda a: sum(list(map(int, a.split(",")))) if a != "." else 0
    )
df
df_filtered_denovo["dad_dp"] = df_filtered_denovo["per_allele_reads_father"].apply(
        lambda a: sum(list(map(int, a.split(",")))) if a != "." else 0
    )
df_filtered_denovo[["per_allele_reads_child","child_coverage","per_allele_reads_father","per_allele_reads_mother","dad_dp","mom_dp"]].head()

# %%
depth_min = 10
df_filtered_denovo = df_filtered_denovo.query(
        f"child_coverage >= {depth_min} and mom_dp >= {depth_min} and dad_dp >= {depth_min}"
    )
df_filtered_denovo.head()
print(len(df_filtered_denovo))
df_filtered_denovo[["allele_origin","denovo_status","per_allele_reads_child","child_coverage","denovo_coverage","per_allele_reads_father","per_allele_reads_mother","dad_dp","mom_dp"]].head()

# %%
parental_overlap_frac_max = 1

df_filtered_denovo[["father_overlap_coverage","mother_overlap_coverage","allele_origin","denovo_status","per_allele_reads_child","child_coverage","denovo_coverage","per_allele_reads_father","per_allele_reads_mother","dad_dp","mom_dp"]].head()

# %%
df_filtered_denovo["dad_ev"] = df_filtered_denovo["father_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))
df_filtered_denovo["mom_ev"] = df_filtered_denovo["mother_overlap_coverage"].apply(lambda o: sum(list(map(int, o.split(",")))))

df_filtered_denovo["dad_ev_frac"] = df_filtered_denovo["dad_ev"] / df_filtered_denovo["dad_dp"]
df_filtered_denovo["mom_ev_frac"] = df_filtered_denovo["mom_ev"] / df_filtered_denovo["mom_dp"]

df_filtered_denovo[["father_overlap_coverage","mother_overlap_coverage","child_coverage","dad_dp","mom_dp","dad_ev","mom_ev","dad_ev_frac","mom_ev_frac"]].head()
# %%
    # calculate the total amount of parental evidence for the de novo AL, aggregated
    # across both parents. calculate the total amount of parental depth at the site.
df_filtered_denovo["parental_ev"] = df_filtered_denovo["dad_ev"] + df_filtered_denovo["mom_ev"]
df_filtered_denovo["parental_dp"] = df_filtered_denovo["dad_dp"] + df_filtered_denovo["mom_dp"]

df_filtered_denovo["parental_ev_frac"] = df_filtered_denovo["parental_ev"] / df_filtered_denovo["parental_dp"]

df_filtered_denovo = df_filtered_denovo[df_filtered_denovo  ["dad_ev_frac"] <= parental_overlap_frac_max]
df_filtered_denovo = df_filtered_denovo[df_filtered_denovo["mom_ev_frac"] <= parental_overlap_frac_max]

df_filtered_denovo = df_filtered_denovo[df_filtered_denovo["parental_ev_frac"] <= parental_overlap_frac_max]

print(len(df_filtered_denovo))
df_filtered_denovo.head()

# %%
from cyvcf2 import VCF
import pandas as pd

def vcf_to_df(filepath):
    records = []
    vcf = VCF(filepath)
    
    for v in vcf:
        # FORMAT fields: use v.format("FIELD") — returns a numpy array (samples × alleles)
        al = v.format("AL")  # shape: (n_samples, n_alleles)
        
        records.append({
            "CHROM":        v.CHROM,
            "POS":          v.POS,
            "END":          v.INFO.get("END"),
            "TRID":         v.INFO.get("TRID"),
            "MOTIFS":       v.INFO.get("MOTIFS"),
            "MOTIF_SIZE":   len(v.INFO.get("MOTIFS", "")),
            "AL_allele1":   int(al[0][0]) if al is not None else None,
            "AL_allele2":   int(al[0][1]) if al is not None and al.shape[1] > 1 else None,
            "AM_allele1":   float(v.format("AM")[0][0]),
            "AM_allele2":   float(v.format("AM")[0][1]),
            "TR_length": int(v.INFO.get("END")) - v.POS 
        })
    
    return pd.DataFrame(records)

df = vcf_to_df("Filter_tr_genotypes/UNMC-000034-01_trgt_filtered_max_sd_ap.vcf")
print(df.head())
print(f"\nTotal STR loci: {len(df):,}")
# %%
merged_df = df.merge(df_filtered_denovo, left_on='TRID', right_on='trid', how='inner')
print(merged_df.head())
print(f"\nTotal de novo STRs with VCF annotations: {len(merged_df):,}")
# %%
merged_df[["CHROM","MOTIFS","AM_allele1","AM_allele2","allele_origin","denovo_status","father_AL","mother_AL","child_AL"]].head()
# %%
import pysam
import pathlib
import sys
import re
merged_df["AM_allele1"].describe()
def GCcontent(sequence):
    """
    Calculate the GC content of a sequence.
    :param sequence: DNA sequence as a string.
    :return: GC content as a float.

    >>> GCcontent('ATTTTGGGGGCCCCCATTTT')
    0.5
    >>> GCcontent('ATTT')
    0.0
    >>> GCcontent('ATTTTGGGGGNNNNNNNNNNNCCCCCATTTT')
    0.5
    >>> GCcontent('GGGGCCC')
    1.0
    >>> GCcontent('') is None
    True
    >>> GCcontent('NNNNNNNNNNN') is None
    True
    """
    if len(sequence) == sequence.count('N'):
        return None
    sequence = sequence.upper()
    return (sequence.count('G') + sequence.count('C')) / (len(sequence) - sequence.count('N'))
# %%
merged_df["GC_CONTENT"] = merged_df["MOTIFS"].apply(GCcontent)
merged_df[["CHROM","MOTIFS","GC_CONTENT","AM_allele1","AM_allele2","allele_origin","denovo_status","father_AL","mother_AL","child_AL"]].head()
merged_df["GC_CONTENT"].describe()
# %%
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np
LABEL_SIZE  = 14
TICK_SIZE   = 12
TITLE_SIZE  = 16
ANNOT_SIZE  = 9
BG          = "white"
GRID_COLOR  = "#e0e0e0"
SPINE_COLOR = "#bbbbbb"

CHR_ORDER = [f"chr{i}" for i in range(1, 24)] + ["chrX", "chrY"]

def sort_chroms(series):
    present = set(series)
    return [c for c in CHR_ORDER if c in present]
# %%
# ══════════════════════════════════════════════════════════════════
# PLOT 2 – Number of TRs per Chromosome
# ══════════════════════════════════════════════════════════════════
canon_df = merged_df[merged_df["CHROM"].isin(CHR_ORDER)]
counts   = canon_df.groupby("CHROM").size()
chroms   = sort_chroms(counts.index)
counts   = counts.reindex(chroms, fill_value=0)

cmap   = plt.cm.get_cmap("tab20b", len(chroms))
colors = [cmap(i) for i in range(len(chroms))]

fig, ax = plt.subplots(figsize=(18, 7))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

xs   = np.arange(len(chroms))
bars = ax.bar(xs, counts.values, width=0.65, color=colors,
              alpha=0.88, linewidth=0, zorder=3)

for bar, cnt in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width() / 2,
            cnt + counts.max() * 0.01,
            f"{cnt:,}", ha="center", va="bottom",
            fontsize=ANNOT_SIZE, color="#333", fontweight="bold")

ax.set_xticks(xs)
ax.set_xticklabels([c.replace("chr", "") for c in chroms],
                   fontsize=TICK_SIZE, color="#222")
ax.set_ylabel("Number of TRs",  fontsize=LABEL_SIZE, color="#222")
ax.set_xlabel("Chromosome",     fontsize=LABEL_SIZE, color="#222")
ax.set_title("Number of Tandem Repeats per Chromosome",
             fontsize=TITLE_SIZE, fontweight="bold", color="#222", pad=12)
ax.tick_params(axis="y", colors="#222", labelsize=TICK_SIZE)
ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
ax.set_ylim(0, counts.max() * 1.12)
ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)

plt.tight_layout()
##plt.savefig("plot_tr_per_chrom.png", dpi=180,
  ##          bbox_inches="tight", facecolor=BG)
plt.show()
# %%
# ══════════════════════════════════════════════════════════════════
# PLOT 3 – Motif Proportion per Chromosome (stacked bar)
# ══════════════════════════════════════════════════════════════════
def classify_motif_size(size):
    if size == 1:    return "1"
    elif size == 2:  return "2"
    elif size == 3:  return "3"
    elif size == 4:  return "4"
    elif size == 5:  return "5"
    elif size == 6:  return "6"
    elif size <= 20: return "7-20"
    else:            return ">20"

CATEGORIES = ["1", "2", "3", "4", "5", "6", "7-20", ">20"]
COLORS = {
    "1":    "#f5e626", "2": "#a0da39", "3":    "#4ac16d",
    "4":    "#1fa187", "5": "#277f8e", "6":    "#365c8d",
    "7-20": "#46327e", ">20":          "#440154",
}

merged_df["MOTIF_CLASS"] = merged_df["MOTIF_SIZE"].apply(classify_motif_size)
plot_df = merged_df[merged_df["CHROM"].isin(CHR_ORDER[:-1])].copy()  # exclude chrM
chroms  = [c for c in CHR_ORDER if c in plot_df["CHROM"].values]

counts_m = (
    plot_df.groupby(["CHROM", "MOTIF_CLASS"])
    .size()
    .unstack(fill_value=0)
    .reindex(index=chroms, columns=CATEGORIES, fill_value=0)
)
proportions = counts_m.div(counts_m.sum(axis=1), axis=0) * 100

fig, ax = plt.subplots(figsize=(18, 7))
fig.patch.set_facecolor(BG)
ax.set_facecolor(BG)

xs      = np.arange(len(chroms))
bar_w   = 0.7
bottoms = np.zeros(len(chroms))

for cat in CATEGORIES:
    vals = proportions[cat].values
    ax.bar(xs, vals, width=bar_w, bottom=bottoms,
           color=COLORS[cat], linewidth=0.2, edgecolor="white")
    bottoms += vals

ax.set_xlim(-0.6, len(chroms) - 0.4)
ax.set_ylim(0, 100)
ax.set_xticks(xs)
ax.set_xticklabels([c.replace("chr", "") for c in chroms],
                   fontsize=TICK_SIZE, color="#222")
ax.set_ylabel("Proportion of TRs (%)", fontsize=LABEL_SIZE, color="#222")
ax.set_xlabel("Chromosome",            fontsize=LABEL_SIZE, color="#222")
ax.set_title("Motif Length Distribution per Chromosome",
             fontsize=TITLE_SIZE, fontweight="bold", color="#222", pad=12)
ax.tick_params(colors="#222", labelsize=TICK_SIZE)
ax.spines[["top", "right"]].set_visible(False)
ax.spines[["left", "bottom"]].set_color(SPINE_COLOR)
ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.7, zorder=0)
ax.set_axisbelow(True)

legend_patches = [mpatches.Patch(color=COLORS[c], label=f"{c} bp")
                  for c in CATEGORIES]
ax.legend(handles=legend_patches, title="Motif length (bp)",
          title_fontsize=TICK_SIZE, fontsize=TICK_SIZE,
          loc="upper right", framealpha=0.9, edgecolor=SPINE_COLOR)

plt.tight_layout()
##plt.savefig("plot_motif_proportion_per_chrom.png", dpi=180,
  ##          bbox_inches="tight", facecolor=BG)
plt.show()

# %%
import seaborn as sns
import matplotlib.pyplot as plt
# Parse allele_origin: extract parent (F/M) ignoring the allele number
merged_df['origin_parent'] = merged_df['allele_origin'].str.extract(r'^([FM])', expand=False).fillna('?')

# Average methylation across both alleles
merged_df['AM_mean'] = merged_df[['AM_allele1', 'AM_allele2']].mean(axis=1)

# %%
sns.violinplot(data=merged_df, x='origin_parent', y='AM_mean', 
               order=['F', 'M', '?'], palette='Set2')
plt.xlabel('Allele Origin (F=Father, M=Mother, ?=Unknown)')
plt.ylabel('Mean Methylation')

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats

df= merged_df.copy()  # use the merged dataframe with annotations
#%%
# --- Classify repeat type based on motif size ---
def classify_repeat(motif_size):
    if motif_size == 1:
        return 'Homopolymer'
    elif 2 <= motif_size <= 6:
        return 'STR'
    else:
        return 'VNTR'

df['repeat_type'] = df['MOTIF_SIZE'].apply(classify_repeat)

df['repeat_type'] = df['MOTIF_SIZE'].apply(classify_repeat)

# --- Parse origin ---
df['origin_parent'] = df['allele_origin'].str.extract(r'^([FM])').fillna('?')[0]
df = df[df['origin_parent'].isin(['F', 'M'])]

# --- Allele length from child_AL (longer allele) ---
df['allele_length'] = df['child_AL'].apply(
    lambda x: max(map(int, str(x).split(','))) if pd.notna(x) else np.nan
)

# --- Count de novo mutations per allele_length bin × origin × repeat_type ---
# Group by repeat type + origin + allele_length and count rows (each row = 1 de novo)
grouped = (df.groupby(['repeat_type', 'origin_parent', 'allele_length'])
             .size()
             .reset_index(name='n_denovo'))

# --- Plot ---
colors = {'F': '#E8A83E', 'M': '#4C8DC4'}
labels = {'F': 'Paternal', 'M': 'Maternal'}
panels = ['Homopolymer', 'STR', 'VNTR']
panel_labels = ['a', 'b', 'c']

fig, axes = plt.subplots(1, 3, figsize=(14, 5))

for ax, panel, letter in zip(axes, panels, panel_labels):
    sub = grouped[grouped['repeat_type'] == panel]

    for origin, grp in sub.groupby('origin_parent'):
        c = colors[origin]
        x = grp['allele_length'].values
        y = grp['n_denovo'].values

        ax.scatter(x, y, color=c, alpha=0.65, s=40,
                   edgecolors='white', linewidths=0.6, zorder=3)

        if len(x) > 2:
            slope, intercept, r, p, se = stats.linregress(x, y)
            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = slope * x_line + intercept

            n = len(x)
            x_mean = x.mean()
            se_line = se * np.sqrt(1/n + (x_line - x_mean)**2 / np.sum((x - x_mean)**2))
            ci = 1.96 * se_line

            ax.plot(x_line, y_line, color=c, linewidth=1.8, zorder=4)
            ax.fill_between(x_line, y_line - ci, y_line + ci,
                            color=c, alpha=0.2, zorder=2)

    ax.set_title(panel, fontsize=12)
    ax.set_xlabel('Allele length', fontsize=10)
    ax.set_ylabel('Number of de novo mutations', fontsize=10)
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(-0.12, 1.05, letter, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top')

# Legend on last panel
patches = [mpatches.Patch(color=colors[o], label=labels[o]) for o in ['M', 'F']]
axes[2].legend(handles=patches, title='Parent-of-origin',
               frameon=False, fontsize=9, title_fontsize=9)

plt.tight_layout()
plt.savefig('denovo_by_repeat_type.png', dpi=150, bbox_inches='tight')
plt.show()

# %%

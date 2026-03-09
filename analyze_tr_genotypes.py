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
            "TR_length": int(v.INFO.get("END")) - v.POS 
        })
    
    return pd.DataFrame(records)

df = vcf_to_df("TR_genotypes/UNMC-000034-01_trgt_genotypes.vcf.gz")
print(df.head())
print(f"\nTotal STR loci: {len(df):,}")
# %%
# ── Overall summary of numeric columns ───────────────────────────
print(df[["TR_length", "AL_allele1", "AL_allele2", "MOTIF_SIZE"]].describe())


# ── Motif size counts ─────────────────────────────────────────────
print("\nMotif Size Distribution:")
print(df["MOTIF_SIZE"].value_counts().sort_index())

#########################
# %%%

from cyvcf2 import VCF
import pandas as pd

def vcf_to_df(filepath):
    records = []
    vcf = VCF(filepath)

    for v in vcf:
        # Skip SNV records — only process STR loci
        if v.INFO.get("VT") != "str":
            continue

        motif    = v.INFO.get("MOTIF", "")
        bed_start = v.INFO.get("BED_START")
        bed_end   = v.INFO.get("BED_END")
        tr_length = (bed_end - bed_start) if bed_start is not None and bed_end is not None else None

        # MC: motif copy number per allele e.g. "3,3"
        mc_raw = None
        try:
            mc_field = v.format("MC")  # numpy array shape (n_samples, 1) — string
            if mc_field is not None:
                mc_raw = mc_field[0][0]  # e.g. "3,3"
        except Exception:
            pass

        mc_allele1, mc_allele2 = None, None
        if mc_raw:
            parts = str(mc_raw).split(",")
            try: mc_allele1 = float(parts[0])
            except: pass
            try: mc_allele2 = float(parts[1])
            except: pass

        # DP: total read depth
        dp = None
        try:
            dp_field = v.format("DP")
            if dp_field is not None:
                dp = int(dp_field[0][0])
        except Exception:
            pass

        # GT
        gt = v.gt_types[0] if v.gt_types is not None else None  # 0=HOM_REF,1=HET,2=UNKNOWN,3=HOM_ALT

        records.append({
            "CHROM":       v.CHROM,
            "POS":         v.POS,
            "ID":          v.ID,
            "BED_START":   bed_start,
            "BED_END":     bed_end,
            "TR_length":   tr_length,
            "MOTIF":       motif,
            "MOTIF_SIZE":  len(motif),
            "REFMC":       v.INFO.get("REFMC"),   # ref motif copy number
            "MC_allele1":  mc_allele1,             # called copy number allele 1
            "MC_allele2":  mc_allele2,             # called copy number allele 2
            "DP":          dp,
            "GT":          gt,
        })

    return pd.DataFrame(records)

df = vcf_to_df("TR_genotypes/UNMC-000034-01_strkit_genotypes.vcf")
print(df.head())
print(f"\nTotal STR loci: {len(df):,}")
# %%
df.shape
# %%
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Shared style settings ─────────────────────────────────────────
LABEL_SIZE  = 14
TICK_SIZE   = 12
TITLE_SIZE  = 16
ANNOT_SIZE  = 9
BG          = "white"
GRID_COLOR  = "#e0e0e0"
SPINE_COLOR = "#bbbbbb"

CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def sort_chroms(series):
    present = set(series)
    return [c for c in CHR_ORDER if c in present]

# ══════════════════════════════════════════════════════════════════
# PLOT 1 – TR Length Distribution (genome-wide)
# ══════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(16, 6))
fig.patch.set_facecolor(BG)
fig.suptitle("TR Length Distribution — Genome Wide",
             fontsize=TITLE_SIZE, fontweight="bold", color="#222")

lengths = df["TR_length"].dropna()
lengths = lengths[lengths > 0]

# Left: full range log scale
ax = axes[0]
ax.set_facecolor(BG)
bins_log = np.logspace(np.log10(lengths.min()), np.log10(lengths.max()), 80)
ax.hist(lengths, bins=bins_log, color="#2196f3", alpha=0.85,
        edgecolor="white", linewidth=0.3)
ax.set_xscale("log")
ax.set_title("Full Range (log scale)", fontsize=TITLE_SIZE - 2, color="#222")
ax.set_xlabel("TR Length (bp)",  fontsize=LABEL_SIZE, color="#222")
ax.set_ylabel("Count",           fontsize=LABEL_SIZE, color="#222")

# Right: zoomed ≤ 200 bp
ax = axes[1]
ax.set_facecolor(BG)
ax.hist(lengths[lengths <= 200], bins=60, color="#43a047", alpha=0.85,
        edgecolor="white", linewidth=0.3)
ax.set_title("Zoomed ≤ 200 bp", fontsize=TITLE_SIZE - 2, color="#222")
ax.set_xlabel("TR Length (bp)", fontsize=LABEL_SIZE, color="#222")
ax.set_ylabel("Count",          fontsize=LABEL_SIZE, color="#222")

for ax in axes:
    ax.tick_params(colors="#222", labelsize=TICK_SIZE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(SPINE_COLOR)
    ax.yaxis.grid(True, color=GRID_COLOR, linewidth=0.7)
    ax.set_axisbelow(True)

plt.tight_layout()
plt.savefig("plot_tr_length_distribution.png", dpi=180,
            bbox_inches="tight", facecolor=BG)
plt.show()

# ══════════════════════════════════════════════════════════════════
# PLOT 2 – Number of TRs per Chromosome
# ══════════════════════════════════════════════════════════════════
canon_df = df[df["CHROM"].isin(CHR_ORDER)]
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
plt.savefig("plot_tr_per_chrom.png", dpi=180,
            bbox_inches="tight", facecolor=BG)
plt.show()

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

df["MOTIF_CLASS"] = df["MOTIF_SIZE"].apply(classify_motif_size)
plot_df = df[df["CHROM"].isin(CHR_ORDER[:-1])].copy()  # exclude chrM
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
plt.savefig("plot_motif_proportion_per_chrom.png", dpi=180,
            bbox_inches="tight", facecolor=BG)
plt.show()

# %%

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
# %%
# ── Overall summary of numeric columns ───────────────────────────
print(df[["TR_length", "AL_allele1", "AL_allele2", "MOTIF_SIZE"]].describe())


# ── Motif size counts ─────────────────────────────────────────────
print("\nMotif Size Distribution:")
print(df["MOTIF_SIZE"].value_counts().sort_index())


# %%
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Chromosome order ──────────────────────────────────────────────
CHR_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

def sort_chroms(series):
    present = set(series)
    return [c for c in CHR_ORDER if c in present]

# ══════════════════════════════════════════════════════════════════
# PLOT 1 – TR Length Distribution (genome-wide)
# ══════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.patch.set_facecolor("#0f0f1a")
fig.suptitle("TR Length Distribution — Genome Wide", fontsize=14,
             color="white", fontweight="bold")

lengths = df["TR_length"].dropna()
lengths = lengths[lengths > 0]

# Left: full range log scale
ax = axes[0]
ax.set_facecolor("#0f0f1a")
bins_log = np.logspace(np.log10(lengths.min()), np.log10(lengths.max()), 80)
ax.hist(lengths, bins=bins_log, color="#4fc3f7", alpha=0.88,
        edgecolor="#0f0f1a", linewidth=0.3)
ax.set_xscale("log")
ax.set_title("Full Range (log scale)", color="white", fontsize=11)
ax.set_xlabel("TR Length (bp)", color="white")
ax.set_ylabel("Count", color="white")

# Right: zoomed ≤ 200 bp
ax = axes[1]
ax.set_facecolor("#0f0f1a")
ax.hist(lengths[lengths <= 200], bins=60, color="#81c784", alpha=0.88,
        edgecolor="#0f0f1a", linewidth=0.3)
ax.set_title("Zoomed ≤ 200 bp", color="white", fontsize=11)
ax.set_xlabel("TR Length (bp)", color="white")
ax.set_ylabel("Count", color="white")

for ax in axes:
    ax.tick_params(colors="white")
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#444")
    ax.yaxis.grid(True, color="#2a2a4a", linewidth=0.5)

plt.tight_layout()
plt.savefig("plot_tr_length_distribution.png", dpi=180,
            bbox_inches="tight", facecolor=fig.get_facecolor())
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

fig, ax = plt.subplots(figsize=(16, 6))
fig.patch.set_facecolor("#0f0f1a")
ax.set_facecolor("#0f0f1a")

xs = np.arange(len(chroms))
bars = ax.bar(xs, counts.values, width=0.65, color=colors,
              alpha=0.88, linewidth=0, zorder=3)

# Value labels on top of each bar
for bar, cnt in zip(bars, counts.values):
    ax.text(bar.get_x() + bar.get_width() / 2,
            cnt + counts.max() * 0.01,
            f"{cnt:,}", ha="center", va="bottom",
            fontsize=6, color="white", fontweight="bold")

ax.set_xticks(xs)
ax.set_xticklabels([c.replace("chr", "") for c in chroms],
                   fontsize=9, color="white")
ax.set_ylabel("Number of TRs", color="white", fontsize=11)
ax.set_title("Number of Tandem Repeats per Chromosome",
             color="white", fontsize=14, fontweight="bold", pad=12)
ax.tick_params(colors="white")
ax.spines[["top", "right", "left", "bottom"]].set_visible(False)
ax.set_ylim(0, counts.max() * 1.12)
ax.yaxis.grid(True, color="#2a2a4a", linewidth=0.5, zorder=0)

plt.tight_layout()
plt.savefig("plot_tr_per_chrom.png", dpi=180,
            bbox_inches="tight", facecolor=fig.get_facecolor())
plt.show()
# %%

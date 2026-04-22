import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── READ FROM XLSX ────────────────────────────────────────────────────────────

XLSX_PATH = "pangenome_data.xlsx"  # ← update path if needed

df_genes  = pd.read_excel(XLSX_PATH, sheet_name="Gene_Lengths")
df_pan    = pd.read_excel(XLSX_PATH, sheet_name="Pangenome_Curves")

core_data = df_genes.loc[df_genes["Gene_type"] == "Core",     "Length_bp"].values
var_data  = df_genes.loc[df_genes["Gene_type"] == "Variable", "Length_bp"].values
n_core    = len(df_genes[df_genes["Gene_type"] == "Core"])
n_var     = len(df_genes[df_genes["Gene_type"] == "Variable"])

# ── PLOT ──────────────────────────────────────────────────────────────────────

CORE_COL = "#4472C4"   # blue
VAR_COL  = "#ED7D31"   # orange

fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))
fig.patch.set_facecolor("white")

# ── Panel A ───────────────────────────────────────────────────────────────────
def jitter(n, width=0.18):
    return np.random.uniform(-width, width, n)

# scatter
ax_a.scatter(0 + jitter(len(core_data)), core_data,
             color=CORE_COL, alpha=0.25, s=2, linewidths=0, zorder=1)
ax_a.scatter(1 + jitter(len(var_data)),  var_data,
             color=VAR_COL,  alpha=0.25, s=2, linewidths=0, zorder=1)

# box-plot helper
def draw_box(ax, x, data, color):
    q1, med, q3 = np.percentile(data, [25, 50, 75])
    iqr = q3 - q1
    lo  = max(data.min(), q1 - 1.5 * iqr)
    hi  = min(data.max(), q3 + 1.5 * iqr)
    mean_v = data.mean()
    w = 0.25
    # filled IQR box
    rect = plt.Rectangle((x - w, q1), 2*w, iqr,
                          facecolor=color, alpha=0.35, edgecolor=color,
                          linewidth=1.5, zorder=3)
    ax.add_patch(rect)
    # median line
    ax.plot([x - w, x + w], [med, med], color=color, linewidth=2, zorder=4)
    # whiskers
    ax.plot([x, x], [lo, q1], color=color, linewidth=1.5, zorder=3)
    ax.plot([x, x], [q3, hi], color=color, linewidth=1.5, zorder=3)
    # caps
    ax.plot([x - w*0.5, x + w*0.5], [lo, lo], color=color, linewidth=1.5, zorder=3)
    ax.plot([x - w*0.5, x + w*0.5], [hi, hi], color=color, linewidth=1.5, zorder=3)
    # mean diamond
    ax.scatter([x], [mean_v], marker="D", color="black", s=30, zorder=5)

draw_box(ax_a, 0, core_data, CORE_COL)
draw_box(ax_a, 1, var_data,  VAR_COL)

# significance bracket
y_top = ax_a.get_ylim()[1] if ax_a.get_ylim()[1] > 1 else 3e4
bh = 6e4
ax_a.annotate("", xy=(1, bh), xytext=(0, bh),
              arrowprops=dict(arrowstyle="-", color="black", lw=1.5))
ax_a.text(0.5, bh * 1.12, "***", ha="center", va="bottom", fontsize=14)

ax_a.set_yscale("log")
ax_a.set_xlim(-0.6, 1.6)
ax_a.set_ylim(90, 2e5)
ax_a.set_xticks([0, 1])
ax_a.set_xticklabels([f"Core genes\n(n={n_core:,})",
                      f"Variable genes\n(n={n_var:,})"], fontsize=11)
ax_a.set_ylabel("Gene length (bp, log scale)", fontsize=11)
ax_a.tick_params(axis="y", labelsize=10)
ax_a.spines[["top","right"]].set_visible(False)
ax_a.text(-0.12, 1.02, "A", transform=ax_a.transAxes,
          fontsize=16, fontweight="bold")

# ── Panel B ───────────────────────────────────────────────────────────────────
gc  = df_pan["Genome_count"].values
cm  = df_pan["Core_mean"].values
cs  = df_pan["Core_std"].values
vm  = df_pan["Variable_mean"].values
vs  = df_pan["Variable_std"].values

ax_b2 = ax_b.twinx()

# Core (left y-axis)
ax_b.plot(gc, cm, color=VAR_COL, linewidth=2, zorder=4)
ax_b.errorbar(gc, cm, yerr=cs, fmt="s", color=VAR_COL,
              markerfacecolor="white", markeredgecolor=VAR_COL,
              markersize=6, capsize=3, linewidth=1, zorder=3)

# Variable (right y-axis)
ax_b2.plot(gc, vm, color=CORE_COL, linewidth=2, zorder=4)
ax_b2.errorbar(gc, vm, yerr=vs, fmt="s", color=CORE_COL,
               markerfacecolor="white", markeredgecolor=CORE_COL,
               markersize=6, capsize=3, linewidth=1, zorder=3)

ax_b.set_xlabel("Number of genomes", fontsize=11)
ax_b.set_ylabel("Core genes", fontsize=11)
ax_b2.set_ylabel("Variable genes", fontsize=11)

ax_b.set_xlim(0.5, 34.5)
ax_b.set_xticks(range(1, 35))
ax_b.tick_params(axis="x", labelsize=7)
ax_b.tick_params(axis="y", labelsize=9)
ax_b2.tick_params(axis="y", labelsize=9)
ax_b.spines[["top"]].set_visible(False)
ax_b2.spines[["top"]].set_visible(False)

legend_handles = [
    mpatches.Patch(color=VAR_COL,  label="Core genes"),
    mpatches.Patch(color=CORE_COL, label="Variable genes"),
]
ax_b.legend(handles=legend_handles, loc="center right", fontsize=10,
            frameon=True, framealpha=0.9)

ax_b.text(-0.12, 1.02, "B", transform=ax_b.transAxes,
          fontsize=16, fontweight="bold")

plt.tight_layout(pad=2.0)
plt.savefig("figure_pangenome.png", dpi=200, bbox_inches="tight")
plt.close()
print("Saved figure_pangenome.png")

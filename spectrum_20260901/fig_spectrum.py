from pathlib import Path
import sys
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
import pandas as pd, numpy as np
import matplotlib.pyplot as plt

setup_plot()
import matplotlib as mpl
mpl.rcParams['axes.facecolor']='white'
mpl.rcParams['figure.facecolor']='white'
mpl.rcParams['axes.grid']=False
res = pd.read_csv("spectrum_mr_results.csv")
order = ["I9_AF","CARDIAC_ARRHYTM","I9_PAROXTAC","I9_OTHARR","I9_AVBLOCK","I9_CONDUCTIO"]
labels = {
 "I9_AF":"Atrial fibrillation and flutter\n(positive control; 63,532 cases)",
 "CARDIAC_ARRHYTM":"Any cardiac arrhythmia (superset;\nincludes AF; 92,926 cases)",
 "I9_PAROXTAC":"Paroxysmal tachycardia (12,878 cases)",
 "I9_OTHARR":"Other arrhythmias (39,265 cases)",
 "I9_AVBLOCK":"Atrioventricular block (7,850 cases)",
 "I9_CONDUCTIO":"Conduction disorders (14,027 cases)",
}
res = res.set_index("endpoint").loc[order].reset_index()

fig, ax = plt.subplots(figsize=(9.2, 5.3))
ypos = np.arange(len(res))[::-1]
for y, (_, r) in zip(ypos, res.iterrows()):
    is_pc = r.endpoint == "I9_AF"
    color = "#B2182B" if is_pc else ("#757575" if r.endpoint=="CARDIAC_ARRHYTM" else "#2166AC")
    ax.plot([r.IVW_lo, r.IVW_hi], [y, y], color=color, lw=2.2, solid_capstyle="round", zorder=2)
    ax.scatter([r.IVW_OR], [y], s=90 if is_pc else 70, color=color,
               marker="s" if is_pc else ("D" if r.endpoint=="CARDIAC_ARRHYTM" else "o"), zorder=3,
               edgecolor="white", linewidth=0.8)
    ax.text(2.62, y, f"OR {r.IVW_OR:.2f} ({r.IVW_lo:.2f}\u2013{r.IVW_hi:.2f}),  P = {r.IVW_P:.3f}".replace("P = 0.", "P = ."),
            va="center", ha="left", fontsize=9, color="#333333")
ax.axvline(1.0, color="#666666", lw=1.0, ls="--", zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([labels[e] for e in res.endpoint], fontsize=9.5)
ax.set_xscale("log")
ax.set_xlim(0.5, 4.2)
ax.set_xticks([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0])
ax.set_xticklabels(["0.5","0.75","1.0","1.5","2.0","3.0","4.0"])
from matplotlib.ticker import NullLocator, NullFormatter
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel("Odds ratio per one cup/day increase in genetically proxied coffee intake (log scale)", fontsize=9.5)
ax.spines[["top","right"]].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_title("Genetic liability to coffee intake across the cardiac arrhythmia spectrum (FinnGen R12, IVW, 40 variants)", fontsize=10.5, pad=10)
fig.savefig("FigureS_Arrhythmia_Spectrum.png", dpi=400, bbox_inches="tight")
fig.savefig("FigureS_Arrhythmia_Spectrum.tiff", dpi=400, bbox_inches="tight", pil_kwargs={"compression":"tiff_lzw"})
print("saved")

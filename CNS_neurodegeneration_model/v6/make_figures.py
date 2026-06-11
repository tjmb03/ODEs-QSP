"""Generate v6 README figures from the actual analysis data."""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import glaucoma_abc_v6 as G
import drug_timing_demo as D

plt.rcParams.update({
    "figure.facecolor": "white", "axes.facecolor": "#fbfbfd",
    "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold",
    "axes.edgecolor": "#cccccc", "axes.grid": True,
    "grid.color": "#e6e6ee", "grid.linewidth": 0.8, "axes.axisbelow": True,
})
CTRL, TRT = "#E0556B", "#2BB6A3"   # disease / therapy
IDC, PARTC = "#3FA34D", "#E0A92B"  # identifiable / partial
DOSES = [90, 150, 210]
IDX = dict(RPE=0, M_mig=3, M1=4, C1q=6, C3a=8, NTF=12, RGC=13)

# ---- run control & treated arms (same 17-state model) ----
tc, Yc = D.run(D.P, [], t_end=365)
tt, Yt = D.run(D.P, DOSES, t_end=365)

def dose_lines(ax):
    for d in DOSES:
        ax.axvline(d, color="#9b8cf0", ls=(0, (4, 3)), lw=1, alpha=0.7, zorder=1)

# === fig1: RGC + RPE survival ===
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for a, idx, ttl in [(ax[0], IDX["RGC"], "RGC survival"), (ax[1], IDX["RPE"], "RPE survival")]:
    a.plot(tc, Yc[idx]*100, color=CTRL, lw=2.4, label="control")
    a.plot(tt, Yt[idx]*100, color=TRT, lw=2.4, label="treated")
    dose_lines(a); a.set_title(ttl); a.set_xlabel("day"); a.set_ylabel("% surviving")
    a.set_xlim(0, 365); a.legend(frameon=False, loc="lower left")
ax[0].annotate(f"+{(Yt[IDX['RGC']][-1]-Yc[IDX['RGC']][-1])*100:.1f} pp",
               (365, Yt[IDX['RGC']][-1]*100), (300, Yt[IDX['RGC']][-1]*100+6),
               color=TRT, fontweight="bold", fontsize=11)
fig.suptitle("v6 disease trajectory — peptide therapy vs natural history", fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig1_endpoints.png", dpi=150, bbox_inches="tight"); plt.close()

# === fig2: complement loop nodes ===
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for a, idx, ttl in [(ax[0], IDX["C1q"], "C1q  (node A)"),
                    (ax[1], IDX["C3a"], "C3a  (node B, loop closure)"),
                    (ax[2], IDX["M1"], "M1 microglia")]:
    a.plot(tc, Yc[idx], color=CTRL, lw=2.4, label="control")
    a.plot(tt, Yt[idx], color=TRT, lw=2.4, label="treated")
    dose_lines(a); a.set_title(ttl); a.set_xlabel("day"); a.set_xlim(0, 365)
    a.legend(frameon=False, loc="upper left")
ax[0].set_ylabel("relative level")
fig.suptitle("Complement feedback loop — node activity, control vs treated", fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig2_loop_nodes.png", dpi=150, bbox_inches="tight"); plt.close()

# === fig3: profile-likelihood curves (the practical-identifiability figure) ===
db = json.load(open("pl_curves_full.json"))
thr = db["thresh"]
order = ["k_C1q_M1","k_C3_cleave","k_C3aR_act","k_damp_M1","k_mig_damp","k_M1_switch",
         "k_des_on","K_anti_M1","k_rpe_cyt","k_rgc_rpe","k_rgc_cyt"]
fig, axes = plt.subplots(3, 4, figsize=(15, 9.5))
for ax_, name in zip(axes.ravel(), order):
    r = db["params"][name]; g = np.array(r["grid"]); d = np.array(r["dchi2"])
    col = IDC if r["status"] == "IDENTIFIABLE" else PARTC
    ax_.plot(g, d, "-o", color=col, ms=4, lw=2)
    ax_.axhline(thr, color="#888", ls="--", lw=1)
    ax_.axvline(r["theta_hat"], color="#555", ls=":", lw=1)
    ax_.text(0.04, 0.92, name, transform=ax_.transAxes, fontweight="bold", fontsize=11, va="top")
    ax_.text(0.96, 0.92, r["status"].replace("IDENTIFIABLE","identif."),
             transform=ax_.transAxes, ha="right", va="top", fontsize=9, color=col, fontweight="bold")
    ax_.set_ylim(-2, min(40, max(8, d.max()*1.1)))
    ax_.set_xlabel("parameter value"); ax_.set_ylabel(r"$\Delta\chi^2$")
axes.ravel()[-1].axis("off")
axes.ravel()[-1].legend(handles=[Patch(color=IDC, label="identifiable (finite CI)"),
                                 Patch(color=PARTC, label="partial (one-sided)"),
                                 plt.Line2D([],[],color="#888",ls="--",label=r"$\chi^2_{0.95,1}=3.84$")],
                        loc="center", frameon=False, fontsize=11)
fig.suptitle("Profile likelihood (pseudo-data) — practical identifiability, v6 braked model",
             fontweight="bold", fontsize=14)
fig.tight_layout(rect=[0,0,1,0.97]); fig.savefig("figures/fig3_profile_likelihood.png", dpi=150, bbox_inches="tight"); plt.close()

# === fig4: singular-value spectrum (structural vs practical rank) ===
si = json.load(open("glaucoma_structural_id_v6.json"))
sv = np.array(si["singular_values_operating"]); sv = sv/sv[0]
fig, ax = plt.subplots(figsize=(9, 4.6))
cols = [IDC if s >= 1e-3 else PARTC for s in sv]
ax.bar(range(1, 12), sv, color=cols, edgecolor="#444", lw=0.6)
ax.set_yscale("log"); ax.set_ylim(1e-5, 2)
ax.axhline(1e-3, color="#888", ls="--", lw=1.2)
ax.text(11.4, 1.2e-3, "practical cutoff (1e-3)", ha="right", fontsize=9, color="#666")
ax.axhline(1e-8, color="#bbb", ls=":", lw=1.2)
ax.set_xticks(range(1, 12)); ax.set_xlabel(r"singular value index"); ax.set_ylabel(r"$\sigma_i/\sigma_{max}$")
ax.set_title("Sensitivity-matrix spectrum at the operating point")
ax.text(0.5, 0.93, "structural rank 11/11 (all real, SIAN-confirmed)   ·   practical rank 9/11",
        transform=ax.transAxes, ha="center", fontsize=10, color="#333")
ax.legend(handles=[Patch(color=IDC, label=r"$\sigma>10^{-3}\sigma_{max}$ (well-determined)"),
                   Patch(color=PARTC, label="weak directions = death-rate params")],
          frameon=False, fontsize=9, loc="lower left")
fig.tight_layout(); fig.savefig("figures/fig4_singular_spectrum.png", dpi=150, bbox_inches="tight"); plt.close()

# === fig5: per-parameter leverage ===
cn = si["column_norms"]; names = sorted(cn, key=cn.get)
partials = {"k_rpe_cyt","k_rgc_rpe","k_rgc_cyt"}
fig, ax = plt.subplots(figsize=(9, 5))
cols = [PARTC if n in partials else IDC for n in names]
ax.barh(range(len(names)), [cn[n] for n in names], color=cols, edgecolor="#444", lw=0.6)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
ax.set_xlabel("output leverage  ‖S column‖"); ax.set_title("Per-parameter leverage on the observables")
ax.legend(handles=[Patch(color=IDC, label="signalling / brake params"),
                   Patch(color=PARTC, label="death-rate params (practically weak)")],
          frameon=False, fontsize=10, loc="lower right")
fig.tight_layout(); fig.savefig("figures/fig5_leverage.png", dpi=150, bbox_inches="tight"); plt.close()

print("figures written:")
import os
for f in sorted(os.listdir("figures")): print("  figures/"+f, f"({os.path.getsize('figures/'+f)//1024} KB)")

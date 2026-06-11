"""Generate validation figures: cross-solver agreement + calibration scorecard."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from scipy.integrate import solve_ivp
import glaucoma_abc_v6 as G
import drug_timing_demo as D

plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"#fbfbfd","font.size":11,
    "axes.titlesize":12,"axes.titleweight":"bold","axes.edgecolor":"#cccccc","axes.grid":True,
    "grid.color":"#e6e6ee","grid.linewidth":0.8,"axes.axisbelow":True})
IDC, PARTC, REF = "#3FA34D", "#E0A92B", "#E0556B"

# ===== fig6: cross-solver agreement (control arm, multiple integrators) =====
P = D.P
C3ss = P['k_C3_base']/P['k_C3_deg']; NTFss = P['k_ntf_base']/P['k_deg_ntf']
y0 = [1,0,P['M_total'],0,0,0,0,C3ss,0,0,0,0,NTFss,1,0,0,0]
grid = np.linspace(0, 365, 366)
methods = {"LSODA":"#2BB6A3","RK45":"#9b8cf0","Radau":"#E0A92B","DOP853":"#E0556B"}
rgc = {}
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios":[2,1]})
for m, c in methods.items():
    sol = solve_ivp(D.odes, (0,365), y0, args=(P,), t_eval=grid, method=m, rtol=1e-7, atol=1e-9)
    rgc[m] = np.interp(grid, sol.t, sol.y[13]*100)
    ax[0].plot(grid, rgc[m], color=c, lw=2.2, alpha=0.85, label=m)
ax[0].set_title("RGC survival under four independent integrators")
ax[0].set_xlabel("day"); ax[0].set_ylabel("% surviving"); ax[0].set_xlim(0,365)
ax[0].legend(frameon=False, loc="lower left")
# max pairwise deviation
stack = np.vstack(list(rgc.values()))
maxdev = (stack.max(0) - stack.min(0))
ax[1].plot(grid, maxdev, color="#444", lw=1.8)
ax[1].set_title("max deviation across solvers")
ax[1].set_xlabel("day"); ax[1].set_ylabel("Δ RGC (pp)"); ax[1].set_xlim(0,365)
ax[1].text(0.5, 0.9, f"peak Δ = {maxdev.max():.2e} pp\n(stiff + explicit agree)",
           transform=ax[1].transAxes, ha="center", va="top", fontsize=10, color="#333")
fig.suptitle("Numerical validation — solver independence (LSODA / RK45 / Radau / DOP853)", fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig6_solver_agreement.png", dpi=150, bbox_inches="tight"); plt.close()

# ===== fig7: calibration scorecard (6 biological targets) =====
s = G.summary_stats(dict(G.PARAM_DEFAULTS))
labels = {"RGC_pct":"RGC survival %","RPE_pct":"RPE survival %","M1_peak":"M1 peak %",
          "ignition_day":"loop ignition (day)","cyt_frac":"cytokine death frac %","trans_frac":"trans-syn death frac %"}
keys = list(labels)
fig, ax = plt.subplots(figsize=(10, 4.8))
for i, k in enumerate(keys):
    lo, hi, _ = G.TARGETS[k]; val = s[k]; norm = (val-lo)/(hi-lo); inside = 0 <= norm <= 1
    ax.add_patch(Rectangle((0, i-0.3), 1, 0.6, color="#d8efdc", ec="#3FA34D", lw=1.2, zorder=1))
    ax.plot(norm, i, "D", color=(IDC if inside else REF), ms=12, zorder=3, mec="#222", mew=0.7)
    ax.text(-0.06, i, labels[k], ha="right", va="center", fontsize=10.5)
    ax.text(1.06, i, f"{val:.1f}  ∈ [{lo:.0f}, {hi:.0f}]  {'✓' if inside else '✗'}",
            ha="left", va="center", fontsize=10, color=(IDC if inside else REF), fontweight="bold")
ax.set_xlim(-0.55, 1.7); ax.set_ylim(-0.7, len(keys)-0.3)
ax.set_xticks([0,1]); ax.set_xticklabels(["target min","target max"], fontsize=9)
ax.set_yticks([]); ax.grid(False)
for sp in ["top","right","left"]: ax.spines[sp].set_visible(False)
ax.set_title("Calibration validation — v6 defaults hit 6 / 6 biological targets", fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig7_calibration_targets.png", dpi=150, bbox_inches="tight"); plt.close()

print("validation figures written: fig6_solver_agreement.png, fig7_calibration_targets.png")
print(f"  solver peak deviation: {maxdev.max():.2e} pp over 365 d")
print(f"  targets inside band: {sum(0<=(s[k]-G.TARGETS[k][0])/(G.TARGETS[k][1]-G.TARGETS[k][0])<=1 for k in keys)}/6")

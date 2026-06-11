"""Stiffness demonstration: explicit fixed-step RK4 vs implicit/adaptive at high
complement-degradation rate, + the Jacobian stiffness ratio."""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import drug_timing_demo as D

plt.rcParams.update({"figure.facecolor":"white","axes.facecolor":"#fbfbfd","font.size":11,
    "axes.titlesize":12,"axes.titleweight":"bold","axes.edgecolor":"#cccccc","axes.grid":True,
    "grid.color":"#e6e6ee","grid.linewidth":0.8,"axes.axisbelow":True})

P = D.P
def odes_unclamped(t, y, p):
    # identical physics to drug_timing_demo.odes but WITHOUT the max(v,0) state clamp
    (RPE,DAMPs,M0,M_mig,M1,M2,C1q,C3,C3a,C5a,Cyt_pro,Cyt_anti,NTF,RGC,A,Cp,R_des) = y
    S=max(0,(p['IOP_target']-p['IOP_normal'])/p['IOP_normal'])
    prot=NTF/(p['EC50_ntf']+NTF)
    rpe_d=p['k_rpe_stress']*S+p['k_rpe_cyt']*Cyt_pro+p['k_rpe_phago']*M_mig
    rgc_d=(p['k_rgc_cyt']*Cyt_pro+p['k_rgc_rpe']*(1-RPE)+p['k_rgc_iop']*S)*(1-prot)
    mig=(p['k_mig_damp']*DAMPs+p['k_mig_C5a']*C5a)*M0
    inhib=1/(1+Cyt_anti/p['K_anti_M1'])
    tlr4=p['k_damp_M1']*DAMPs*M_mig*inhib
    c3ar=p['k_C3aR_act']*C3a*(1-R_des)*M_mig*inhib
    M1sw=p['k_M1_switch']*M1
    C3cl=p['k_C3_cleave']*C1q*C3
    return [-rpe_d*RPE,
            p['k_damp_rpe']*rpe_d*RPE+p['k_damp_rgc']*rgc_d*RGC-p['k_damp_clear']*DAMPs,
            -mig+p['k_deact_M1']*M1+p['k_res_M2']*M2+p['k_return']*M_mig,
            mig-c3ar-tlr4-p['k_return']*M_mig, tlr4+c3ar-M1sw-p['k_deact_M1']*M1, M1sw-p['k_res_M2']*M2,
            p['k_C1q_M1']*M1-p['k_C1q_deg']*C1q, p['k_C3_base']+p['k_C3_rpe']*rpe_d*RPE-C3cl-p['k_C3_deg']*C3,
            p['k_C3a_frac']*C3cl-p['k_C3a_deg']*C3a, p['k_C5a_frac']*C3cl-p['k_C5a_deg']*C5a,
            p['k_M1_cyt']*M1-p['k_deg_pro']*Cyt_pro-p['k_inhib']*Cyt_anti*Cyt_pro,
            p['k_M2_cyt']*M2-p['k_deg_anti']*Cyt_anti, p['k_ntf_base']+p['k_M2_ntf']*M2-p['k_deg_ntf']*NTF,
            -rgc_d*RGC, -p['k_abs']*A, p['k_abs']*A-p['k_el_pep']*Cp,
            p['k_des_on']*C3a*(1-R_des)-p['k_des_off']*R_des]

def rk4_fixed(f, y0, p, t_end, dt, clamp):
    n = int(t_end/dt); y = np.array(y0, float); ts=[0.0]; ys=[y.copy()]
    for i in range(n):
        t = i*dt
        if clamp: y = np.maximum(y, 0.0)
        k1=np.array(f(t,y,p)); k2=np.array(f(t,y+dt/2*k1,p))
        k3=np.array(f(t,y+dt/2*k2,p)); k4=np.array(f(t,y+dt*k3,p))
        y = y + dt/6*(k1+2*k2+2*k3+k4)
        ts.append((i+1)*dt); ys.append(y.copy())
    return np.array(ts), np.array(ys)

# 1) get an active-loop state (run control to day 200)
tc, Yc = D.run(P, [], t_end=200)
y_mid = Yc[:, -1].copy()

# 2) make complement fast/stiff
pstiff = dict(P); pstiff['k_C3a_deg'] = 500.0; pstiff['k_C5a_deg'] = 500.0
IDX_C3A = 8
tau = 1/pstiff['k_C3a_deg']                       # 0.002 d fast timescale
dt_lim = 2.78*tau                                 # RK4 abs-stability limit ~0.0056 d
T = 0.3

ref = solve_ivp(odes_unclamped, (0,T), y_mid, args=(pstiff,), method='LSODA',
                rtol=1e-9, atol=1e-12, dense_output=True)
tg = np.linspace(0, T, 600); c3a_ref = ref.sol(tg)[IDX_C3A]
t_big, Y_big = rk4_fixed(odes_unclamped, y_mid, pstiff, T, 0.02, clamp=False)   # dt > limit
t_bigc, Y_bigc = rk4_fixed(D.odes,       y_mid, pstiff, T, 0.02, clamp=True)    # clamped masks it
t_sm, Y_sm  = rk4_fixed(odes_unclamped, y_mid, pstiff, T, 0.002, clamp=False)   # dt < limit

# 3) stiffness ratio via finite-diff Jacobian eigenvalues
def stiff_ratio(p, y):
    n=len(y); J=np.zeros((n,n)); f0=np.array(odes_unclamped(0,y,p)); h=1e-7
    for j in range(n):
        yp=y.copy(); yp[j]+=h; J[:,j]=(np.array(odes_unclamped(0,yp,p))-f0)/h
    ev=np.linalg.eigvals(J); re=np.abs(ev.real); re=re[re>1e-9]
    return re.max()/re.min(), re
r_def,_ = stiff_ratio(P, y_mid)
r_stf,_ = stiff_ratio(pstiff, y_mid)

# ---- plot ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"width_ratios":[2,1]})
eps = 1e-30
peak_big = np.nanmax(np.abs(Y_big[:,IDX_C3A]))
ax[0].semilogy(tg, np.abs(c3a_ref)+eps, color="#2BB6A3", lw=3, label="LSODA (adaptive/implicit) — reference", zorder=2)
ax[0].semilogy(t_sm, np.abs(Y_sm[:,IDX_C3A])+eps, "--", color="#3FA34D", lw=1.9, label="explicit RK4, dt=0.002 (< stability limit) ✓")
ax[0].semilogy(t_big, np.abs(Y_big[:,IDX_C3A])+eps, "-o", color="#E0556B", ms=3, lw=1.9, label="explicit RK4, dt=0.02 (> limit) — UNSTABLE", zorder=3)
ax[0].semilogy(t_bigc, np.abs(Y_bigc[:,IDX_C3A])+eps, ":", color="#E0A92B", lw=2, label="same dt, state-clamped — bounded but biased")
ax[0].set_title(f"|C3a| at high degradation (k_C3a_deg=500, τ≈{tau:.3f} d)\nstiff mode: explicit RK4 unstable for dt > {dt_lim:.4f} d")
ax[0].set_xlabel("day (short window)"); ax[0].set_ylabel("|C3a|  (log scale)"); ax[0].set_xlim(0, T)
ax[0].set_ylim(1e-4, 1e12)
ax[0].annotate(f"diverges → {peak_big:.0e}", xy=(0.12, 1e9), fontsize=9.5, color="#E0556B", fontweight="bold")
ax[0].legend(frameon=False, fontsize=8.5, loc="lower right")

ax[1].bar(["default\nrates","stiff\n(k_deg=500)"], [r_def, r_stf], color=["#2BB6A3","#E0556B"], edgecolor="#444")
ax[1].set_yscale("log"); ax[1].set_title("Jacobian stiffness ratio\n|λ|max / |λ|min")
ax[1].set_ylabel("stiffness ratio (log)")
for i,v in enumerate([r_def, r_stf]):
    ax[1].text(i, v*1.3, f"{v:.0f}×", ha="center", fontweight="bold", fontsize=11)
fig.suptitle("Numerical stiffness — why fixed-step explicit RK4 needs help (motivates ROS2 / QSSA)", fontweight="bold")
fig.tight_layout(); fig.savefig("figures/fig8_stiffness.png", dpi=150, bbox_inches="tight"); plt.close()
print(f"fig8_stiffness.png written")
print(f"  default stiffness ratio: {r_def:.0f}x | stiff: {r_stf:.0f}x")
print(f"  RK4 abs-stability limit at k_deg=500: dt < {dt_lim:.4f} d; dt=0.02 is {0.02/dt_lim:.1f}x over -> unstable")
print(f"  explicit dt=0.02 peak |C3a|: {np.nanmax(np.abs(Y_big[:,IDX_C3A])):.2e}  vs reference max {c3a_ref.max():.2e}")

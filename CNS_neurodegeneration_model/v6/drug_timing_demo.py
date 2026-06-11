"""Demonstrate: dosing aligned to the cascade ignition window restores efficacy.
Adds drug PK (A_eye, C_pep) + 4 PD Hill terms to the v6 model, matching the R app."""
import numpy as np
from scipy.integrate import solve_ivp

# Report's actual params
P = dict(IOP_normal=15, IOP_target=21, M_total=1.0,
  k_rpe_stress=0.003, k_rpe_cyt=0.006, k_rpe_phago=0.006,
  k_damp_rpe=0.40, k_damp_rgc=0.20, k_damp_clear=0.80,
  k_mig_damp=1.85, k_mig_C5a=0.80, k_return=0.08,
  k_damp_M1=0.50, k_C3aR_act=2.28, k_M1_switch=0.176, k_deact_M1=0.05, k_res_M2=0.12,
  k_C1q_M1=1.36, k_C1q_deg=0.40, k_C3_base=0.30, k_C3_cleave=1.20, k_C3_deg=0.30, k_C3_rpe=0.20,
  k_C3a_frac=0.60, k_C3a_deg=0.50, k_C5a_frac=0.30, k_C5a_deg=0.40,
  k_M1_cyt=1.0, k_deg_pro=0.35, k_inhib=0.25, k_M2_cyt=0.70, k_deg_anti=0.28,
  k_ntf_base=0.05, k_M2_ntf=1.80, k_deg_ntf=0.28,
  k_rgc_cyt=0.018, k_rgc_rpe=0.005, k_rgc_iop=0.002, EC50_ntf=0.40,
  k_des_on=5.0, k_des_off=0.111, K_anti_M1=0.595,
  k_abs=1.0, k_el_pep=0.10,   # ~6.9-day half-life (realistic intravitreal Fab)
  Emax_C1qblock=6.0, Emax_C3aRblock=5.0, Emax_switch=7.0, Emax_migration=5.0,
  EC50_pep=0.80, gamma_pep=2.0)

def hill(C,Emax,EC50,g):
    if C<=0: return 0.0
    return Emax*C**g/(EC50**g+C**g)

def odes(t,y,p):
    (RPE,DAMPs,M0,M_mig,M1,M2,C1q,C3,C3a,C5a,Cyt_pro,Cyt_anti,NTF,RGC,A,Cp,R_des)=[max(v,0) for v in y]
    RPE=min(RPE,1);RGC=min(RGC,1);R_des=min(R_des,1)
    S=max(0,(p['IOP_target']-p['IOP_normal'])/p['IOP_normal'])
    PC1q=hill(Cp,p['Emax_C1qblock'],p['EC50_pep'],p['gamma_pep'])
    PC3aR=hill(Cp,p['Emax_C3aRblock'],p['EC50_pep'],p['gamma_pep'])
    Psw=hill(Cp,p['Emax_switch'],p['EC50_pep'],p['gamma_pep'])
    Pmig=hill(Cp,p['Emax_migration'],p['EC50_pep'],p['gamma_pep'])
    rpe_d=p['k_rpe_stress']*S+p['k_rpe_cyt']*Cyt_pro+p['k_rpe_phago']*M_mig
    dRPE=-rpe_d*RPE
    prot=NTF/(p['EC50_ntf']+NTF)
    rgc_d=(p['k_rgc_cyt']*Cyt_pro+p['k_rgc_rpe']*(1-RPE)+p['k_rgc_iop']*S)*(1-prot)
    dRGC=-rgc_d*RGC
    dDAMPs=p['k_damp_rpe']*rpe_d*RPE+p['k_damp_rgc']*rgc_d*RGC-p['k_damp_clear']*DAMPs
    mig=(p['k_mig_damp']*DAMPs+p['k_mig_C5a']*C5a)*(1-Pmig)*M0
    Mret=p['k_return']*M_mig
    dM0=-mig+p['k_deact_M1']*M1+p['k_res_M2']*M2+Mret
    inhib=1/(1+Cyt_anti/p['K_anti_M1'])
    tlr4=p['k_damp_M1']*DAMPs*M_mig*inhib
    c3ar=p['k_C3aR_act']*C3a*(1-R_des)*M_mig*(1-PC3aR)*inhib
    M1sw=p['k_M1_switch']*(1+Psw)*M1
    dM_mig=mig-c3ar-tlr4-Mret
    dM1=tlr4+c3ar-M1sw-p['k_deact_M1']*M1
    dM2=M1sw-p['k_res_M2']*M2
    dR_des=p['k_des_on']*C3a*(1-R_des)-p['k_des_off']*R_des
    dC1q=p['k_C1q_M1']*M1*(1-PC1q)-p['k_C1q_deg']*C1q
    C3cl=p['k_C3_cleave']*C1q*C3
    dC3=p['k_C3_base']+p['k_C3_rpe']*rpe_d*RPE-C3cl-p['k_C3_deg']*C3
    dC3a=p['k_C3a_frac']*C3cl-p['k_C3a_deg']*C3a
    dC5a=p['k_C5a_frac']*C3cl-p['k_C5a_deg']*C5a
    dCp_pro=p['k_M1_cyt']*M1-p['k_deg_pro']*Cyt_pro-p['k_inhib']*Cyt_anti*Cyt_pro
    dCa=p['k_M2_cyt']*M2-p['k_deg_anti']*Cyt_anti
    dNTF=p['k_ntf_base']+p['k_M2_ntf']*M2-p['k_deg_ntf']*NTF
    dA=-p['k_abs']*A
    dCp=p['k_abs']*A-p['k_el_pep']*Cp
    return [dRPE,dDAMPs,dM0,dM_mig,dM1,dM2,dC1q,dC3,dC3a,dC5a,dCp_pro,dCa,dNTF,dRGC,dA,dCp,dR_des]

def run(p,doses,t_end=365):
    C3ss=p['k_C3_base']/p['k_C3_deg']; NTFss=p['k_ntf_base']/p['k_deg_ntf']
    y=[1,0,p['M_total'],0,0,0,0,C3ss,0,0,0,0,NTFss,1, (p['dose_amount'] if 0 in doses else 0),0,0]
    bp=sorted([d for d in doses if d>0])+[t_end]
    ts=[]; ys=[]; tc=0
    for tb in bp:
        sol=solve_ivp(odes,(tc,tb),y,args=(p,),t_eval=np.linspace(tc,tb,max(50,int((tb-tc)*4))),method='LSODA',rtol=1e-7,atol=1e-9)
        ts.append(sol.t); ys.append(sol.y)
        y=sol.y[:,-1].tolist()
        if tb in doses: y[14]+=p['dose_amount']  # A_eye dose
        tc=tb
    t=np.concatenate(ts); Y=np.concatenate(ys,axis=1)
    return t,Y

def finalRGC(t,Y): return Y[13,-1]*100

P['dose_amount']=5.0
# Decompose the fix: PK half-life vs dosing schedule. Base P uses the SHIPPING
# default k_el_pep=0.10 (~6.9-day half-life); we override to 1.5 to reproduce the
# report's fast-clearing peptide where the original problem appeared.
FAST = dict(P, k_el_pep=1.5)   # report's PK (t1/2 = 0.46 d) — drug gone by ~day 35
REAL = dict(P)                 # shipping PK (t1/2 = 6.9 d)

t,Yc   = run(REAL, [])                          # control (no drug)
t,Yfe  = run(FAST, [0,7,14,21,28])              # report's setup: fast PK + early dosing
t,Yre  = run(REAL, [0,7,14,21,28])              # fix A: realistic PK only (still early dosing)
t,Yfw  = run(FAST, [90,120,150,180,210,240])    # fix B: window dosing only (still fast PK)
t,Yrw  = run(REAL, [90,150,210])                # fix C (shipping default): realistic PK + window

base=finalRGC(t,Yc)
def row(lbl,Y): print(f"  {lbl:<46} {finalRGC(t,Y):>5.1f}%   (Δ {finalRGC(t,Y)-base:+.1f}pp)")
print("RGC survival at day 365  (control = %.1f%%):\n" % base)
row("Report: fast PK (t½=0.5d) + early dosing d0-28", Yfe)
row("Fix A — realistic PK (t½=6.9d), same early dosing", Yre)
row("Fix B — window dosing d90-240, fast PK", Yfw)
row("Fix C — realistic PK + window d90/150/210 [DEFAULT]", Yrw)
print("\nTakeaway: the drug is dosed before the cascade ignites (~day 156), so a")
print("fast-clearing peptide is gone before it matters (+small effect). Either a")
print("longer half-life OR dosing in the disease window restores efficacy; the")
print("shipping default does both (3 bi-monthly doses, realistic intravitreal PK).")

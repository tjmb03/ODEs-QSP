# glaucoma_v6_SIAN.jl — GLOBAL structural identifiability, v6 braked model
# Constants baked in as rationals; ONLY the 11 calibration parameters are symbolic:
#   k_C1q_M1 k_C3_cleave k_C3aR_act k_damp_M1 k_rpe_cyt k_rgc_rpe k_rgc_cyt
#   k_mig_damp k_M1_switch k_des_on K_anti_M1
# Outputs y1..y5 = RGC, RPE, M1, C1q, C3a (control arm; drug/PK states omitted).
#
# SETUP (once):
#   curl -fsSL https://install.julialang.org | sh      # installs Julia via juliaup
#   julia
#   julia> using Pkg; Pkg.add("StructuralIdentifiability")
#
# RUN:
#   julia> include("glaucoma_v6_SIAN.jl")
#   julia> assess_local_identifiability(ode)           # FAST — start here (Sedoglavic)
#   julia> assess_identifiability(ode; prob_threshold=0.99)   # GLOBAL — may be slow for 15 states
#   julia> find_identifiable_functions(ode)            # reveals identifiable param COMBINATIONS
#
# TRACTABILITY: 15 states × 11 params rational system. Local ID should complete in
# seconds–minutes. GLOBAL ID (Groebner/differential elimination) can be very slow or
# blow up memory at this size — if it stalls, reduce: keep only the complement-loop
# states (M1,Mmig,C1q,C3,C3a,Rdes) + outputs, or fix more parameters as numeric.

using StructuralIdentifiability

ode = @ODEmodel(
    RPE'(t)   = -((3//1000)*(2//5) + k_rpe_cyt*Cytpro(t) + (3//500)*Mmig(t))*RPE(t),
    DAMPs'(t) =  (2//5)*((3//1000)*(2//5) + k_rpe_cyt*Cytpro(t) + (3//500)*Mmig(t))*RPE(t)
               + (1//5)*(k_rgc_cyt*Cytpro(t) + k_rgc_rpe*(1 - RPE(t)) + (1//500)*(2//5))*RGC(t)
               - (4//5)*DAMPs(t),
    M0'(t)    = -(k_mig_damp*DAMPs(t) + (4//5)*C5a(t))*M0(t)
               + (1//20)*M1(t) + (3//25)*M2(t) + (2//25)*Mmig(t),
    Mmig'(t)  =  (k_mig_damp*DAMPs(t) + (4//5)*C5a(t))*M0(t)
               - k_C3aR_act*C3a(t)*(1 - Rdes(t))*Mmig(t)*(1/(1 + Cytanti(t)/K_anti_M1))
               - k_damp_M1*DAMPs(t)*Mmig(t)*(1/(1 + Cytanti(t)/K_anti_M1))
               - (2//25)*Mmig(t),
    M1'(t)    =  k_damp_M1*DAMPs(t)*Mmig(t)*(1/(1 + Cytanti(t)/K_anti_M1))
               + k_C3aR_act*C3a(t)*(1 - Rdes(t))*Mmig(t)*(1/(1 + Cytanti(t)/K_anti_M1))
               - k_M1_switch*M1(t) - (1//20)*M1(t),
    M2'(t)    =  k_M1_switch*M1(t) - (3//25)*M2(t),
    Rdes'(t)  =  k_des_on*C3a(t)*(1 - Rdes(t)) - (111//1000)*Rdes(t),
    C1q'(t)   =  k_C1q_M1*M1(t) - (2//5)*C1q(t),
    C3'(t)    =  (3//10) + (1//5)*((3//1000)*(2//5) + k_rpe_cyt*Cytpro(t) + (3//500)*Mmig(t))*RPE(t)
               - k_C3_cleave*C1q(t)*C3(t) - (3//10)*C3(t),
    C3a'(t)   =  (3//5)*k_C3_cleave*C1q(t)*C3(t) - (1//2)*C3a(t),
    C5a'(t)   =  (3//10)*k_C3_cleave*C1q(t)*C3(t) - (2//5)*C5a(t),
    Cytpro'(t)=  1*M1(t) - (7//20)*Cytpro(t) - (1//4)*Cytanti(t)*Cytpro(t),
    Cytanti'(t)= (7//10)*M2(t) - (7//25)*Cytanti(t),
    NTF'(t)   =  (1//20) + (9//5)*M2(t) - (7//25)*NTF(t),
    RGC'(t)   = -(k_rgc_cyt*Cytpro(t) + k_rgc_rpe*(1 - RPE(t)) + (1//500)*(2//5))
                 *(1 - NTF(t)/((2//5) + NTF(t)))*RGC(t),
    y1(t) = RGC(t),
    y2(t) = RPE(t),
    y3(t) = M1(t),
    y4(t) = C1q(t),
    y5(t) = C3a(t)
)

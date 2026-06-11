# glaucoma_v6_SIAN_4out.jl - TEST: does v6 structural identifiability survive WITHOUT
# the C3a output? Same full 15-state braked model, but observing only the v5.5 output
# set {RGC, RPE, M1, C1q}. C3a remains a hidden (unobserved) intermediate state.
#
# WHY: the v5.5 SIAN/DAISY analysis reported some parameters structurally
# non-identifiable from {RPE, RGC, M1, C1q}. The v6 model (with the two saturating
# brakes) is structurally identifiable from {RGC,RPE,M1,C1q,C3a}. This file isolates
# the cause: if assess_local_identifiability STILL returns all-true here (4 outputs,
# no C3a), then the v6 identifiability comes from the BRAKES / extra observable
# dynamics, NOT from adding the C3a output -> confirms the README "vs v5.5" claim.
# (Numerical sensitivity-rank proxy already predicts 11/11 even without C3a.)
#
# RUN:
#   julia> include("glaucoma_v6_SIAN_4out.jl")
#   julia> assess_local_identifiability(ode)      # expect: all 11 params + 15 states => true
# Check the summary prints  Outputs: y1, y2, y3, y4  (no y5).

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
    y4(t) = C1q(t)
)

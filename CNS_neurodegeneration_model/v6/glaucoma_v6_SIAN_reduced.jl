# glaucoma_v6_SIAN_reduced.jl - GLOBAL structural identifiability of the DEATH-RATE
# subsystem (the practically-weak directions in the full model).
#
# The full 15-state model is LOCALLY identifiable (assess_local_identifiability = all
# true) but its GLOBAL IO-equation step OOMs. This reduced model isolates the three
# death-rate parameters and runs GLOBAL in seconds, giving the discrete-alias check.
#
# Reduction: keep states RPE, RGC. Treat the three drivers that carry the death terms
# - cytokine (Cytpro), migratory microglia (Mmig), neurotrophic factor (NTF) - as
# KNOWN, independently-measured INPUTS. This is exactly the experimental condition that
# breaks the operating-point collinearity (Cytpro measured independently of RPE loss).
# Outputs y1=RGC, y2=RPE. Parameters assessed: k_rpe_cyt, k_rgc_cyt, k_rgc_rpe.
# All other coefficients baked in as rationals (S, k_rpe_stress, k_rpe_phago,
# k_rgc_iop, EC50_ntf).
#
# RUN:
#   julia> include("glaucoma_v6_SIAN_reduced.jl")
#   julia> assess_identifiability(ode; prob_threshold=0.99)   # GLOBAL - seconds here
#   julia> find_identifiable_functions(ode)
# Check the printed summary shows  Inputs: Cytpro, Mmig, NTF.
#
# EXPECTED: with the cytokine known as an input, Cytpro and (1-RPE) are independent
# signals, so k_rgc_cyt (on Cytpro) and k_rgc_rpe (on 1-RPE) separate -> all three
# :globally identifiable. That is the algebraic statement of "the degeneracy is
# practical: it disappears once cytokine is measured/perturbed independently of RPE."

using StructuralIdentifiability

ode = @ODEmodel(
    RPE'(t) = -(3//2500 + k_rpe_cyt*Cytpro(t) + 3//500*Mmig(t))*RPE(t),
    RGC'(t) = -(k_rgc_cyt*Cytpro(t) + k_rgc_rpe*(1 - RPE(t)) + 1//1250)
               *(1 - NTF(t)/(2//5 + NTF(t)))*RGC(t),
    y1(t) = RGC(t),
    y2(t) = RPE(t)
)

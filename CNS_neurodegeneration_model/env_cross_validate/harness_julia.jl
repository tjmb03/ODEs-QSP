#!/usr/bin/env julia
# harness_julia.jl — Julia engine for cross-validation, SEGMENT-RESTART dosing.
#
# Doses are applied by restarting the integration at each dose boundary
# (u[A_eye] += dose at the segment start, so a dose day's output is POST-dose),
# exactly matching the Python and R engines. This deliberately avoids
# PresetTimeCallback: with its default save_positions=(true,true) it saves BOTH
# pre- and post-dose at each dose time (duplicate rows on the daily grid), and
# the saveat-vs-affect! ordering of the recorded value is version-sensitive —
# the same dose-day ambiguity that bit deSolve. Segment-restart is unambiguous.
#
# Solver is FBDF() (multistep BDF, analog of lsoda/ode15s), independent of the
# other engines. Do NOT use LSODA.jl (it wraps the same Fortran as R's deSolve).
#
# Usage:  julia harness_julia.jl <control|treated> <out.csv>
# Deps:   julia> ] add OrdinaryDiffEq

using OrdinaryDiffEq

const SNAMES = ["RPE","DAMPs","M0","M_mig","M1","M2","C1q","C3","C3a","C5a",
                "Cyt_pro","Cyt_anti","NTF","RGC","A_eye","C_pep","R_des"]

const P = (
    IOP_normal=15.0, IOP_target=21.0, M_total=1.0,
    k_rpe_stress=0.003, k_rpe_cyt=0.006, k_rpe_phago=0.006,
    k_damp_rpe=0.40, k_damp_rgc=0.20, k_damp_clear=0.80,
    k_mig_damp=1.850, k_mig_C5a=0.80, k_return=0.08, k_damp_M1=0.50,
    k_C3aR_act=2.275, k_M1_switch=0.176, k_deact_M1=0.05, k_res_M2=0.12,
    k_C1q_M1=1.359, k_C1q_deg=0.40,
    k_C3_base=0.30, k_C3_cleave=1.20, k_C3_deg=0.30, k_C3_rpe=0.20,
    k_C3a_frac=0.60, k_C3a_deg=0.50, k_C5a_frac=0.30, k_C5a_deg=0.40,
    k_M1_cyt=1.00, k_deg_pro=0.35, k_inhib=0.25, k_M2_cyt=0.70, k_deg_anti=0.28,
    k_ntf_base=0.05, k_M2_ntf=1.80, k_deg_ntf=0.28,
    k_rgc_cyt=0.0179, k_rgc_rpe=0.005, k_rgc_iop=0.002, EC50_ntf=0.40,
    k_des_on=5.000, k_des_off=0.111, K_anti_M1=0.595,
    k_abs=1.00, k_el_pep=0.10,
    Emax_C1qblock=6.0, Emax_C3aRblock=5.0, Emax_switch=7.0, Emax_migration=5.0,
    EC50_pep=0.80, gamma_pep=2.00,
)

hill(C, Emax, EC50, g) = (Cc = max(C, 0.0); Emax * Cc^g / (EC50^g + Cc^g))

function glaucoma!(du, u, p, t)
    RPE = clamp(u[1], 0.0, 1.0); DAMPs = max(u[2], 0.0); M0 = max(u[3], 0.0)
    M_mig = max(u[4], 0.0); M1 = max(u[5], 0.0); M2 = max(u[6], 0.0)
    C1q = max(u[7], 0.0); C3 = max(u[8], 0.0); C3a = max(u[9], 0.0); C5a = max(u[10], 0.0)
    Cp = max(u[11], 0.0); Ca = max(u[12], 0.0); NTF = max(u[13], 0.0); RGC = clamp(u[14], 0.0, 1.0)
    A = max(u[15], 0.0); Cp2 = max(u[16], 0.0); R_des = clamp(u[17], 0.0, 1.0)
    Stress = max(0.0, (p.IOP_target - p.IOP_normal) / p.IOP_normal)
    Pb = hill(Cp2, p.Emax_C1qblock, p.EC50_pep, p.gamma_pep)
    Pc = hill(Cp2, p.Emax_C3aRblock, p.EC50_pep, p.gamma_pep)
    Ps = hill(Cp2, p.Emax_switch, p.EC50_pep, p.gamma_pep)
    Pm = hill(Cp2, p.Emax_migration, p.EC50_pep, p.gamma_pep)
    rpe_d = p.k_rpe_stress * Stress + p.k_rpe_cyt * Cp + p.k_rpe_phago * M_mig
    prot = NTF / (p.EC50_ntf + NTF)
    rgc_d = (p.k_rgc_cyt * Cp + p.k_rgc_rpe * (1 - RPE) + p.k_rgc_iop * Stress) * (1 - prot)
    mig = (p.k_mig_damp * DAMPs + p.k_mig_C5a * C5a) * (1 - Pm) * M0
    Mret = p.k_return * M_mig
    inhib = 1.0 / (1.0 + Ca / p.K_anti_M1)
    tlr4 = p.k_damp_M1 * DAMPs * M_mig * inhib
    c3ar = p.k_C3aR_act * C3a * (1 - R_des) * M_mig * (1 - Pc) * inhib
    M1sw = p.k_M1_switch * (1 + Ps) * M1
    C3cl = p.k_C3_cleave * C1q * C3
    du[1]  = -rpe_d * RPE
    du[2]  = p.k_damp_rpe * rpe_d * RPE + p.k_damp_rgc * rgc_d * RGC - p.k_damp_clear * DAMPs
    du[3]  = -mig + p.k_deact_M1 * M1 + p.k_res_M2 * M2 + Mret
    du[4]  = mig - c3ar - tlr4 - Mret
    du[5]  = tlr4 + c3ar - M1sw - p.k_deact_M1 * M1
    du[6]  = M1sw - p.k_res_M2 * M2
    du[7]  = p.k_C1q_M1 * M1 * (1 - Pb) - p.k_C1q_deg * C1q
    du[8]  = p.k_C3_base + p.k_C3_rpe * rpe_d * RPE - C3cl - p.k_C3_deg * C3
    du[9]  = p.k_C3a_frac * C3cl - p.k_C3a_deg * C3a
    du[10] = p.k_C5a_frac * C3cl - p.k_C5a_deg * C5a
    du[11] = p.k_M1_cyt * M1 - p.k_deg_pro * Cp - p.k_inhib * Ca * Cp
    du[12] = p.k_M2_cyt * M2 - p.k_deg_anti * Ca
    du[13] = p.k_ntf_base + p.k_M2_ntf * M2 - p.k_deg_ntf * NTF
    du[14] = -rgc_d * RGC
    du[15] = -p.k_abs * A
    du[16] = p.k_abs * A - p.k_el_pep * Cp2
    du[17] = p.k_des_on * C3a * (1 - R_des) - p.k_des_off * R_des
    return nothing
end

const DOSE_TIMES  = [0.0, 90.0, 150.0, 210.0]   # CANONICAL — identical across engines
const DOSE_AMOUNT = 5.0

function run_arm(arm)
    C3_ss  = P.k_C3_base  / P.k_C3_deg
    NTF_ss = P.k_ntf_base / P.k_deg_ntf
    u0 = zeros(17)
    u0[1] = 1.0; u0[3] = P.M_total; u0[8] = C3_ss; u0[13] = NTF_ss; u0[14] = 1.0
    # A_eye (u[15]) starts at 0; the t=0 loading dose is applied in the loop.

    ds       = arm == "treated" ? DOSE_TIMES : Float64[]
    bounds   = sort(unique(vcat(0.0, ds, 365.0)))
    all_days = collect(0.0:1.0:365.0)

    rec_t = Float64[]
    rec_u = Vector{Vector{Float64}}()
    u = copy(u0)
    for i in 1:(length(bounds) - 1)
        t0 = bounds[i]; t1 = bounds[i + 1]
        if t0 in ds
            u[15] += DOSE_AMOUNT                    # dose at boundary -> day is POST-dose
        end
        seg = filter(d -> t0 <= d < t1, all_days)   # output days [t0, t1)
        if i == length(bounds) - 1
            push!(seg, t1)                          # final segment: include endpoint
        end
        sa  = sort(unique(vcat(seg, t1)))           # integrate THROUGH t1
        sol = solve(ODEProblem(glaucoma!, u, (t0, t1), P), FBDF();
                    saveat = sa, reltol = 1e-8, abstol = 1e-10)
        segset = Set(seg)
        for (j, t) in enumerate(sol.t)
            if t in segset
                push!(rec_t, t); push!(rec_u, copy(sol.u[j]))
            end
        end
        u = copy(sol.u[end])                        # carry state at t1 to next dose
    end
    return rec_t, rec_u
end

arm = length(ARGS) >= 1 ? ARGS[1] : "control"
out = length(ARGS) >= 2 ? ARGS[2] : "julia_engine.csv"
rec_t, rec_u = run_arm(arm)

open(out, "w") do io
    println(io, "time," * join(SNAMES, ","))
    for (t, uu) in zip(rec_t, rec_u)
        println(io, string(t) * "," * join(string.(uu), ","))
    end
end
println("Julia/FBDF (segment-restart dosing) wrote $(length(rec_t)) rows x " *
        "$(length(SNAMES)) states ($arm arm) to $out")

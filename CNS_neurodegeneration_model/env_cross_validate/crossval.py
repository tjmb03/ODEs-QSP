"""
crossval.py — cross-engine validation across R / Python / MATLAB.

The HTML/JS dashboard is the interactive front-end, NOT a validation target.
Scientific cross-validation compares three independent stiff-capable engines:

    R       deSolve::lsoda          (harness_r.R, sources glaucoma_qsp_v6.R)
    Python  scipy LSODA             (this module)
    MATLAB  SimBiology ode15s       (harness_matlab.m)

Each engine integrates the SAME canonical spec (params, ICs, dose schedule,
time grid) and writes a trajectory CSV with NAMED state columns. This
comparator matches columns BY NAME — so the R-vs-Python C1q/C3 index-ordering
difference is irrelevant — interpolates every engine onto the canonical daily
grid, and reports pairwise agreement plus any single outlier engine.

Protocol: cross-validate the CONTROL arm first (isolates the vector field from
the dosing-schedule ambiguity), then the TREATED arm with a schedule replicated
identically in all three harnesses.
"""
from __future__ import annotations
from dataclasses import dataclass
import csv, itertools
import numpy as np
import validation_tools as vt

CANON_GRID = np.arange(0, 366, 1.0)          # daily, day 0..365
CANON_DOSES = (0, 90, 150, 210)              # explicit; all engines must match


# ── Python engine: integrate the canonical spec, write a named-column CSV ───
def export_python_trajectory(path, arm="control", method="radau",
                             overrides=None, dose_times=CANON_DOSES):
    """Python engine. Default solver is Radau (implicit RK) — deliberately NOT
    LSODA, because scipy's LSODA shares the ODEPACK Fortran lineage with R's
    deSolve::lsoda; Radau keeps Python independent of the R engine."""
    sys = vt.build_glaucoma_system(treated=(arm == "treated"),
                                   dose_times=dose_times)
    if overrides:
        sys = sys.with_params(**overrides)
    t, Y = vt.integrate(sys, method=method, dt=0.05)
    cols = {s: np.interp(CANON_GRID, t, Y[:, sys.idx(s)]) for s in sys.state_names}
    _write_csv(path, CANON_GRID, cols)
    return path


def _write_csv(path, grid, cols):
    names = list(cols.keys())
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["time"] + names)
        for i, t in enumerate(grid):
            w.writerow([t] + [cols[n][i] for n in names])


def load_trajectory(path):
    """Load a CSV (time + named state cols) and interpolate onto CANON_GRID,
    so engines that emit different output grids still align."""
    with open(path) as f:
        r = csv.DictReader(f)
        names = [n for n in r.fieldnames if n.lower() != "time"]
        rows = list(r)
    t = np.array([float(x["time"]) for x in rows])
    order = np.argsort(t); t = t[order]
    return {n: np.interp(CANON_GRID, t,
                         np.array([float(x[n]) for x in rows])[order]) for n in names}


@dataclass(frozen=True)
class CrossValResult:
    passed: bool
    engines: list
    readouts: list
    pairwise_readout_dev: dict     # "A|B" -> {readout: max_abs_dev}
    worst: dict                    # readout -> (pair, dev, ratio_to_tol)
    outlier: str | None
    summary: str


def cross_validate(engine_csvs, readouts=("RGC", "RPE", "M1"),
                   atol=1e-5, rtol=5e-3) -> CrossValResult:
    """engine_csvs: {engine_name: csv_path}. Verdict: every pairwise readout
    deviation within atol + rtol*scale (max over the whole trajectory). The
    0.5% floor reflects what independent stiff adaptive solvers (lsoda / LSODA /
    ode15s) realistically agree to; tighten to ~1e-3 once you confirm the three
    real engines beat it (they should — RK4 is the loose one, and it's excluded)."""
    traj = {e: load_trajectory(p) for e, p in engine_csvs.items()}
    engines = list(traj)
    shared = set.intersection(*[set(t) for t in traj.values()])
    miss = [r for r in readouts if r not in shared]
    if miss:
        raise ValueError(f"readout(s) {miss} not present in all engines")
    pairwise, worst = {}, {}
    disagree = {e: 0 for e in engines}
    all_ok = True
    for a, b in itertools.combinations(engines, 2):
        key = f"{a}|{b}"; pairwise[key] = {}
        for s in readouts:
            dev = float(np.max(np.abs(traj[a][s] - traj[b][s])))
            thr = atol + rtol * float(np.max(np.abs(traj[a][s])))
            pairwise[key][s] = dev
            if dev > thr:
                all_ok = False; disagree[a] += 1; disagree[b] += 1
            if s not in worst or dev > worst[s][1]:
                worst[s] = (key, dev, dev / max(thr, 1e-12))
    outlier = None
    if not all_ok:
        m = max(disagree.values())
        tops = [e for e, c in disagree.items() if c == m]
        outlier = tops[0] if len(tops) == 1 else None
    msg = (f"{len(engines)} engines agree on {list(readouts)} "
           f"(atol={atol:g}, rtol={rtol:g})" if all_ok else
           "readout disagreement" +
           (f"; likely outlier: {outlier}" if outlier else "; no single outlier"))
    return CrossValResult(all_ok, engines, list(readouts), pairwise, worst,
                          outlier, msg + " -> " + ("PASS" if all_ok else "FAIL"))


def conservation_from_csv(path, pool=("M0", "M_mig", "M1", "M2"),
                          expected=1.0, tol=1e-6):
    """Engine-agnostic mass-balance check: the microglia pool sum must stay at
    `expected` along the trajectory. Works on any engine's CSV output."""
    traj = load_trajectory(path)
    missing = [c for c in pool if c not in traj]
    if missing:
        raise ValueError(f"pool columns {missing} absent in {path}")
    s = sum(traj[c] for c in pool)
    drift = float(np.max(np.abs(s - expected)))
    return {"passed": drift <= tol, "max_drift": drift, "tol": tol,
            "pool": list(pool), "expected": expected}


# ── Self-test: prove the comparator catches a cross-engine discrepancy ──────
if __name__ == "__main__":
    import tempfile, os
    d = tempfile.mkdtemp()
    good_lsoda = export_python_trajectory(os.path.join(d, "py_lsoda.csv"),
                                          arm="control", method="ref")
    good_rk4 = export_python_trajectory(os.path.join(d, "py_rk4.csv"),
                                        arm="control", method="rk4")
    # simulate a transcription error in one "engine": k_C3_cleave off by 5%
    buggy = export_python_trajectory(os.path.join(d, "buggy_engine.csv"),
                                     arm="control", method="ref",
                                     overrides={"k_C3_cleave": 1.20 * 1.05})

    print("== two correct engines (LSODA vs RK4) ==")
    r1 = cross_validate({"Python_LSODA": good_lsoda, "AltEngine_RK4": good_rk4})
    print("  ", r1.summary)
    for s, (pair, dev, ratio) in r1.worst.items():
        print(f"     {s}: max dev {dev:.2e} ({ratio:.2f}×tol) [{pair}]")

    print("== one engine with a 5% k_C3_cleave transcription error ==")
    r2 = cross_validate({"Python_LSODA": good_lsoda, "AltEngine_RK4": good_rk4,
                         "BuggyEngine": buggy})
    print("  ", r2.summary)
    for s, (pair, dev, ratio) in r2.worst.items():
        print(f"     {s}: max dev {dev:.2e} ({ratio:.2f}×tol) [{pair}]")

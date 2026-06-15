#!/usr/bin/env python3
"""
run_crossval_gate.py — deterministic cross-engine validation GATE (no LLM).

Compares R / Python / MATLAB / Julia trajectory CSVs for the glaucoma v6 model.
Built for CI: exits 0 iff every available engine pair agrees on the readout
states (control AND treated) and every engine conserves the microglia pool.

Modes:
  (default)        compare CSVs already in --csv-dir, named {engine}_{arm}.csv
                   (CI: each engine runs as its own job and uploads its CSVs)
  --run-engines    first run the harnesses for whichever engine binaries are on
                   PATH (Rscript / julia / matlab / python), both arms, into
                   --csv-dir, then compare (Docker / local all-in-one)

Exit codes:  0 = pass   1 = validation failure   2 = setup / <2 engines
"""
from __future__ import annotations
import argparse, glob, os, shutil, subprocess, sys
import crossval as cv

R_SOURCE = os.environ.get("XV_R_SOURCE", "glaucoma_model.R")
ARMS = ("control", "treated")

CMD = {
    "python": lambda arm, out: [sys.executable, "harness_python.py", arm, out],
    "r":      lambda arm, out: ["Rscript", "harness_r.R", arm, out, R_SOURCE],
    "julia":  lambda arm, out: ["julia", "harness_julia.jl", arm, out],
    "matlab": lambda arm, out: ["matlab", "-batch", f"harness_matlab('{arm}','{out}')"],
}
BIN = {"python": sys.executable, "r": "Rscript", "julia": "julia", "matlab": "matlab"}


def available():
    return [e for e, b in BIN.items() if e == "python" or shutil.which(b)]


def run_engines(csv_dir):
    os.makedirs(csv_dir, exist_ok=True)
    for eng in available():
        for arm in ARMS:
            out = os.path.join(csv_dir, f"{eng}_{arm}.csv")
            print(f"[run] {eng}/{arm}")
            r = subprocess.run(CMD[eng](arm, out), capture_output=True,
                               text=True, timeout=1800)
            print(f"  {'ok -> ' + out if r.returncode == 0 else '!! FAILED'}")
            if r.returncode != 0:
                print(r.stderr[-1000:])


def discover(csv_dir):
    found = {}
    for path in glob.glob(os.path.join(csv_dir, "*_*.csv")):
        eng, arm = os.path.basename(path)[:-4].rsplit("_", 1)
        if arm in ARMS:
            found.setdefault(eng, {})[arm] = path
    return found


def gate(csv_dir):
    found = discover(csv_dir)
    out, ok = [], True

    out.append("## Per-engine conservation (microglia pool)")
    for eng in sorted(found):
        for arm in ARMS:
            p = found[eng].get(arm)
            if not p:
                continue
            try:
                c = cv.conservation_from_csv(p)
                ok &= c["passed"]
                out.append(f"- `{eng}/{arm}` drift {c['max_drift']:.2e} -> "
                           f"{'PASS' if c['passed'] else 'FAIL'}")
            except Exception as e:
                ok = False
                out.append(f"- `{eng}/{arm}` ERROR: {e}")

    for arm in ARMS:
        here = {e: found[e][arm] for e in found if arm in found[e]}
        out.append(f"\n## Cross-validation — {arm} ({len(here)} engines: "
                   f"{', '.join(sorted(here)) or 'none'})")
        if len(here) < 2:
            ok = False
            out.append(f"- INSUFFICIENT engines for {arm} (need >=2)")
            continue
        res = cv.cross_validate(here)
        ok &= res.passed
        out.append(f"- {res.summary}")
        if res.outlier:
            out.append(f"- likely outlier engine: **{res.outlier}**")
        for s, (pair, dev, ratio) in res.worst.items():
            out.append(f"  - {s}: max dev {dev:.2e} ({ratio:.2f}×tol) [{pair}]")

    report = "\n".join(out)
    print(report)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a") as f:
            f.write(f"# Cross-engine validation: {'PASS ✅' if ok else 'FAIL ❌'}\n\n"
                    f"{report}\n")
    return ok, len(found)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv-dir", default="csvs")
    ap.add_argument("--run-engines", action="store_true")
    args = ap.parse_args()

    if args.run_engines:
        run_engines(args.csv_dir)

    ok, n_engines = gate(args.csv_dir)
    if n_engines < 2:
        print(f"\n::error::only {n_engines} engine(s) found in {args.csv_dir}; "
              f"need >=2 to cross-validate", file=sys.stderr)
        sys.exit(2)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

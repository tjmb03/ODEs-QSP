#!/usr/bin/env python3
"""diagnose_treated.py — localize where two engine CSVs diverge on an arm.
Usage: python diagnose_treated.py python_treated.csv r_treated.csv"""
import sys, numpy as np, crossval as cv
a = cv.load_trajectory(sys.argv[1]); b = cv.load_trajectory(sys.argv[2])
shared = sorted(set(a) & set(b))
print(f"{'state':9s} {'max|dev|':>10s} {'@day':>5s} {'1st day>1e-7':>13s}")
for s in shared:
    d = np.abs(a[s] - b[s]); k = int(np.argmax(d))
    fd = next((i for i in range(len(d)) if d[i] > 1e-7), -1)
    print(f"{s:9s} {d.max():10.2e} {k:5d} {fd:13d}")

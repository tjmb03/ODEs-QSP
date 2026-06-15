#!/usr/bin/env python3
"""harness_python.py — Python (scipy Radau) engine CLI for cross-validation.

Symmetry with harness_r.R / harness_matlab.m / harness_julia.jl, so every
engine is invoked the same way: <harness> <control|treated> <out.csv>.
"""
import sys
import crossval as cv


def main():
    arm = sys.argv[1] if len(sys.argv) > 1 else "control"
    out = sys.argv[2] if len(sys.argv) > 2 else "python_engine.csv"
    cv.export_python_trajectory(out, arm=arm)        # default solver: Radau
    print(f"Python/Radau wrote trajectory ({arm} arm) to {out}")


if __name__ == "__main__":
    main()

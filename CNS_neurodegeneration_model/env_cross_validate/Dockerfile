# Reproducible cross-engine validation image: R + Python + Julia.
#
# MATLAB is deliberately NOT in this image — it is license-gated and cannot be
# freely containerized; run the MATLAB engine on a self-hosted / traditionally
# licensed runner. The gate degrades gracefully: with R + Python + Julia present
# it cross-validates those three (>=2 engines required).
#
# Build:  docker build -t glaucoma-xval .
# Run:    docker run --rm glaucoma-xval         # runs all 3 engines + the gate
#
# Note: the Julia layer precompiles OrdinaryDiffEq, which is slow and makes the
# image large (expect a multi-minute build). That cost is paid once at build so
# the first solve at run time is fast.

FROM julia:1.10 AS julia

FROM python:3.12-slim

# Julia from the official image
COPY --from=julia /usr/local/julia /usr/local/julia
ENV PATH="/usr/local/julia/bin:${PATH}"

# R + deSolve from Debian's precompiled binary (r-cran-desolve).
# Avoids a from-source CRAN compile, which silently no-ops in this slim image
# (no gfortran): install.packages() only warns on failure, so Rscript exits 0
# and the missing package would surface only at runtime. apt fails loudly.
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base r-cran-desolve ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work

# Python deps
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# R deps: deSolve is already installed (r-cran-desolve, above). Verify it loads
# so an absent package fails the BUILD here rather than at run time.
RUN Rscript -e 'suppressMessages(library(deSolve)); cat("deSolve", as.character(packageVersion("deSolve")), "OK\n")'

# Julia deps (precompile at build time)
COPY harness_julia.jl .
RUN julia -e 'using Pkg; Pkg.add(["OrdinaryDiffEq","DiffEqCallbacks"]); using OrdinaryDiffEq, DiffEqCallbacks'

# Application: kernel + comparator + harnesses + gate
COPY validation_tools.py crossval.py harness_python.py harness_r.R \
     glaucoma_model.R run_crossval_gate.py ./

# Run every engine present (R + Python + Julia here), both arms, then the gate.
CMD ["python", "run_crossval_gate.py", "--run-engines", "--csv-dir", "/work/csvs"]

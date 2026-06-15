"""
validation_agent.py — Agentic CROSS-ENGINE validation: R / Python / MATLAB / Julia.

Scope: the HTML/JS dashboard is the interactive front-end and is NOT validated
here. The agent drives four independent stiff-capable engines on the same
canonical spec, checks each engine's mass balance, then cross-validates their
trajectories pairwise and flags any outlier.

  R       deSolve::lsoda          via  harness_r.R   (sources glaucoma_qsp_v6.R)
  Python  scipy LSODA             via  crossval.export_python_trajectory
  MATLAB  SimBiology ode15s       via  harness_matlab.m
  Julia   OrdinaryDiffEq FBDF     via  harness_julia.jl

Trust contract (unchanged): the model PLANS and INTERPRETS; every number and
verdict comes from the deterministic tools. The SDK bundles the Claude Code CLI
and bills against your plan — no API key; just an authenticated `claude` login.

Run:  python validation_agent.py
"""
from __future__ import annotations
import asyncio, json, os, shutil, subprocess
from claude_agent_sdk import (tool, create_sdk_mcp_server, ClaudeAgentOptions,
                              query, AssistantMessage, TextBlock)
import crossval as cv

# --- configure for your machine -------------------------------------------
WORKDIR     = os.environ.get("XV_WORKDIR", "/tmp/xval")
HARNESS_DIR = os.environ.get("XV_HARNESS_DIR", os.getcwd())
R_SOURCE    = os.environ.get("XV_R_SOURCE", "glaucoma_qsp_v6.R")
os.makedirs(WORKDIR, exist_ok=True)
# ---------------------------------------------------------------------------


def _text(obj) -> dict:
    return {"content": [{"type": "text", "text": json.dumps(obj, default=float)}]}


async def _run(cmd, cwd=None):
    return await asyncio.to_thread(
        lambda: subprocess.run(cmd, capture_output=True, text=True,
                               timeout=900, cwd=cwd))


@tool("run_python_engine",
      "Run the Python (scipy LSODA) engine on the given arm; writes a "
      "trajectory CSV and returns its path. arm in {control,treated}.",
      {"arm": str})
async def run_python_engine(args):
    arm = args.get("arm", "control")
    out = os.path.join(WORKDIR, f"python_{arm}.csv")
    await asyncio.to_thread(cv.export_python_trajectory, out, arm)
    return _text({"engine": "Python", "arm": arm, "csv": out})


@tool("run_r_engine",
      "Run the R (deSolve::lsoda) engine on the given arm via harness_r.R; "
      "writes a trajectory CSV. Returns an error if Rscript is unavailable.",
      {"arm": str})
async def run_r_engine(args):
    arm = args.get("arm", "control")
    if shutil.which("Rscript") is None:
        return _text({"error": "Rscript not on PATH — R engine unavailable"})
    out = os.path.join(WORKDIR, f"r_{arm}.csv")
    res = await _run(["Rscript", os.path.join(HARNESS_DIR, "harness_r.R"),
                      arm, out, R_SOURCE])
    if res.returncode != 0:
        return _text({"error": "R engine failed", "stderr": res.stderr[-800:]})
    return _text({"engine": "R", "arm": arm, "csv": out, "log": res.stdout.strip()})


@tool("run_matlab_engine",
      "Run the MATLAB (SimBiology ode15s) engine on the given arm via "
      "harness_matlab.m; writes a trajectory CSV. Returns an error if matlab "
      "is unavailable.",
      {"arm": str})
async def run_matlab_engine(args):
    arm = args.get("arm", "control")
    if shutil.which("matlab") is None:
        return _text({"error": "matlab not on PATH — MATLAB engine unavailable"})
    out = os.path.join(WORKDIR, f"matlab_{arm}.csv")
    cmd = f"harness_matlab('{arm}','{out}')"
    res = await _run(["matlab", "-batch", cmd], cwd=HARNESS_DIR)
    if res.returncode != 0:
        return _text({"error": "MATLAB engine failed", "stderr": res.stderr[-800:]})
    return _text({"engine": "MATLAB", "arm": arm, "csv": out, "log": res.stdout.strip()})


@tool("run_julia_engine",
      "Run the Julia (OrdinaryDiffEq, FBDF stiff solver) engine on the given "
      "arm via harness_julia.jl; writes a trajectory CSV. Returns an error if "
      "julia is unavailable.",
      {"arm": str})
async def run_julia_engine(args):
    arm = args.get("arm", "control")
    if shutil.which("julia") is None:
        return _text({"error": "julia not on PATH — Julia engine unavailable"})
    out = os.path.join(WORKDIR, f"julia_{arm}.csv")
    res = await _run(["julia", os.path.join(HARNESS_DIR, "harness_julia.jl"),
                      arm, out])
    if res.returncode != 0:
        return _text({"error": "Julia engine failed", "stderr": res.stderr[-800:]})
    return _text({"engine": "Julia", "arm": arm, "csv": out, "log": res.stdout.strip()})


@tool("check_engine_conservation",
      "Engine-agnostic mass-balance check on one engine's CSV: the microglia "
      "pool (M0+M_mig+M1+M2) must stay constant. Catches transcription errors "
      "in any single engine independently of the others.",
      {"csv_path": str})
async def check_engine_conservation(args):
    r = await asyncio.to_thread(cv.conservation_from_csv, args["csv_path"])
    return _text({"check": "conservation", **r})


@tool("cross_validate_engines",
      "Compare engine trajectory CSVs pairwise on the clinical readout states "
      "(RGC/RPE/M1) over the whole trajectory; flags any single outlier engine. "
      "Pass engines as a JSON object mapping engine_name -> csv_path.",
      {"engines_json": str, "arm": str})
async def cross_validate_engines(args):
    engines = json.loads(args["engines_json"])
    r = await asyncio.to_thread(cv.cross_validate, engines)
    return _text({"check": "cross_validation", "arm": args.get("arm"),
                  "passed": r.passed, "summary": r.summary, "outlier": r.outlier,
                  "worst": {k: {"pair": v[0], "dev": v[1], "ratio": v[2]}
                            for k, v in r.worst.items()},
                  "pairwise_readout_dev": r.pairwise_readout_dev})


VALIDATORS = create_sdk_mcp_server(
    name="xval", version="1.0.0",
    tools=[run_python_engine, run_r_engine, run_matlab_engine, run_julia_engine,
           check_engine_conservation, cross_validate_engines])

SYSTEM_PROMPT = """You are a cross-engine ODE validation agent for a QSP portfolio.
Goal: confirm the glaucoma v6 model is implemented consistently across the R
(lsoda), Python (LSODA), MATLAB (SimBiology ode15s), and Julia (OrdinaryDiffEq
FBDF) engines. The HTML/JS dashboard is the front-end and is OUT OF SCOPE —
never treat it as an engine.

HARD RULES:
- Base EVERY conclusion on tool outputs. Never compute or assert a number/verdict.
- Conservation is per-engine integrity (transcription/mass-balance). It is
  necessary but NOT sufficient; it does not establish that two engines agree.

PROTOCOL:
1. CONTROL arm first (isolates the vector field from the dosing ambiguity):
   run each available engine (python, r, matlab, julia) on arm=control. If an
   engine reports unavailable, note it and proceed with the engines you have
   (need >=2 to cross-validate).
2. Run check_engine_conservation on each CSV produced.
3. cross_validate_engines with a JSON map of the control CSVs. If it flags an
   outlier, that engine has a discrepancy — report which engine and which
   readout/pair drives it.
4. Only if control passes, repeat for the TREATED arm (all engines share the
   canonical {0,90,150,210} schedule; a treated-only disagreement implicates
   the dose handling, not the equations).

Finish with: which engines were compared, per-engine conservation, the pairwise
readout agreement (control and treated), any outlier, and a bottom-line verdict
on whether the implementations are equivalent."""

ALLOWED = [f"mcp__xval__{t}" for t in
           ("run_python_engine", "run_r_engine", "run_matlab_engine",
            "run_julia_engine", "check_engine_conservation",
            "cross_validate_engines")]


async def main():
    opts = ClaudeAgentOptions(
        mcp_servers={"xval": VALIDATORS}, allowed_tools=ALLOWED,
        system_prompt=SYSTEM_PROMPT, model="opus", max_turns=30)
    prompt = ("Cross-validate the glaucoma v6 model across the R, Python, "
              "MATLAB, and Julia engines (control arm first, then treated). "
              "Give me the report.")
    async for msg in query(prompt=prompt, options=opts):
        if isinstance(msg, AssistantMessage):
            for b in msg.content:
                if isinstance(b, TextBlock):
                    print(b.text, end="", flush=True)
    print()


if __name__ == "__main__":
    asyncio.run(main())

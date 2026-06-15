function harness_matlab(arm, outfile)
% harness_matlab — MATLAB/SimBiology engine for cross-validation.
%
% Reuses build_glaucoma_simbiology() to construct the model, applies the
% CANONICAL dose schedule, integrates with ode15s on a daily grid, and writes
% a named-column trajectory CSV matching crossval.py. Column order differs from
% R/Python (C1q before C3) but the comparator matches BY NAME, so it aligns.
%
% Usage:
%   harness_matlab('control','matlab_control.csv')
%   harness_matlab('treated','matlab_treated.csv')
%
% Requires build_glaucoma_simbiology.m on the MATLAB path (same folder).

if nargin < 1, arm = 'control'; end
if nargin < 2, outfile = 'matlab_engine.csv'; end

m = build_glaucoma_simbiology();          % nargout>=1 -> returns model, no auto-sim

doseTimes  = [0 90 150 210];              % CANONICAL — identical across all engines
doseAmount = 5.0;

cs = getconfigset(m);
cs.SolverType                      = 'ode15s';
cs.StopTime                        = 365;
cs.SolverOptions.AbsoluteTolerance = 1e-10;
cs.SolverOptions.RelativeTolerance = 1e-8;
cs.SolverOptions.MaxStep           = 1;
cs.SolverOptions.OutputTimes       = (0:1:365)';

if strcmp(arm, 'treated')
    dose            = sbiodose('pep_dose', 'schedule');
    dose.TargetName = 'A_eye';
    dose.Time       = doseTimes;
    dose.Amount     = repmat(doseAmount, size(doseTimes));
    sd = sbiosimulate(m, cs, [], dose);
else
    sd = sbiosimulate(m, cs, [], []);     % control: no dosing
end

names = sd.DataNames(:)';                 % logged species names
T = array2table([sd.Time, sd.Data], 'VariableNames', [{'time'}, names]);
writetable(T, outfile);
fprintf('MATLAB/ode15s wrote %d rows x %d states (%s arm) to %s\n', ...
        height(T), numel(names), arm, outfile);
end

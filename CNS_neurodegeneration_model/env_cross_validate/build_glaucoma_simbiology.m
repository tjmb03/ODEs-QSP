function m = build_glaucoma_simbiology()
% BUILD_GLAUCOMA_SIMBIOLOGY  Construct the glaucoma QSP v6 model in SimBiology
% as a 17-state rate-rule system, faithful to glaucoma_qsp_v6.R / the HTML
% dashboard. Run with no output to simulate + plot, or capture m to inspect.
%
%   m = build_glaucoma_simbiology();
%
% Requires: MATLAB + SimBiology toolbox. Tested against the documented v6
% operating point (RGC ~64% treated at the dashboard dose schedule).

% ========================================================================
% STEP 0.  Model + compile options
%   Unit conversion / dimensional analysis OFF so every quantity is a bare
%   number, matching the dimensionless R reference exactly.
% ========================================================================
m  = sbiomodel('GlaucomaQSP_v6');
cs = getconfigset(m);
cs.CompileOptions.UnitConversion      = false;
cs.CompileOptions.DimensionalAnalysis = false;

% ========================================================================
% STEP 1.  Compartment (capacity 1 -> amount == concentration numerically)
% ========================================================================
comp = addcompartment(m, 'retina', 1);

% ========================================================================
% STEP 2.  Constant parameters (the kinetic rate constants)
%   Stored as a struct, then added in a loop. Edit values here only.
% ========================================================================
pv = struct( ...
  'IOP_normal',15, 'IOP_target',21, 'M_total',1.0, ...
  'k_rpe_stress',0.003, 'k_rpe_cyt',0.006, 'k_rpe_phago',0.006, ...
  'k_damp_rpe',0.40, 'k_damp_rgc',0.20, 'k_damp_clear',0.80, ...
  'k_mig_damp',1.850, 'k_mig_C5a',0.80, 'k_return',0.08, ...
  'k_damp_M1',0.50, ...
  'k_C3aR_act',2.275, ...
  'k_M1_switch',0.176, 'k_deact_M1',0.05, 'k_res_M2',0.12, ...
  'k_C1q_M1',1.359, 'k_C1q_deg',0.40, ...
  'k_C3_base',0.30, 'k_C3_cleave',1.20, 'k_C3_deg',0.30, 'k_C3_rpe',0.20, ...
  'k_C3a_frac',0.60, 'k_C3a_deg',0.50, ...
  'k_C5a_frac',0.30, 'k_C5a_deg',0.40, ...
  'k_M1_cyt',1.00, 'k_deg_pro',0.35, 'k_inhib',0.25, ...
  'k_M2_cyt',0.70, 'k_deg_anti',0.28, ...
  'k_ntf_base',0.05, 'k_M2_ntf',1.80, 'k_deg_ntf',0.28, ...
  'k_rgc_cyt',0.0179, 'k_rgc_rpe',0.005, 'k_rgc_iop',0.002, ...
  'EC50_ntf',0.40, ...
  'k_des_on',5.000, 'k_des_off',0.111, 'K_anti_M1',0.595, ...
  'k_abs',1.00, 'k_el_pep',0.10, ...
  'Emax_C1qblock',6.0, 'Emax_C3aRblock',5.0, ...
  'Emax_switch',7.0, 'Emax_migration',5.0, ...
  'EC50_pep',0.80, 'gamma_pep',2.00);

fn = fieldnames(pv);
for i = 1:numel(fn)
    addparameter(m, fn{i}, pv.(fn{i}));   % ConstantValue = true by default
end

% Steady-state initial conditions (match R: C3 and NTF start at synth/deg SS)
C3_ss  = pv.k_C3_base  / pv.k_C3_deg;     % = 1.0
NTF_ss = pv.k_ntf_base / pv.k_deg_ntf;    % = 0.17857

% ========================================================================
% STEP 3.  Species (17 states) with initial amounts.
%   A_eye starts at 0; the t=0 loading dose is delivered by the dose object
%   so the schedule is explicit and unambiguous.
% ========================================================================
spec = {
  'RPE',      1
  'DAMPs',    0
  'M0',       pv.M_total
  'M_mig',    0
  'M1',       0
  'M2',       0
  'C1q',      0
  'C3',       C3_ss
  'C3a',      0
  'C5a',      0
  'Cyt_pro',  0
  'Cyt_anti', 0
  'NTF',      NTF_ss
  'RGC',      1
  'A_eye',    0
  'C_pep',    0
  'R_des',    0 };
for i = 1:size(spec,1)
    addspecies(comp, spec{i,1}, spec{i,2});
end

% ========================================================================
% STEP 4.  Intermediate variables as non-constant parameters + repeated
%          assignment rules. SimBiology resolves their dependency order.
%   Hill blocks use power(); max(0,.) on Stress via piecewise (portable
%   across versions). Hill's defensive max(C,0) is dropped: C_pep stays >=0.
% ========================================================================
inter = {
 'Stress'      , 'piecewise((IOP_target-IOP_normal)/IOP_normal, IOP_target > IOP_normal, 0)'
 'Pb'          , 'Emax_C1qblock *power(C_pep,gamma_pep)/(power(EC50_pep,gamma_pep)+power(C_pep,gamma_pep))'
 'Pc'          , 'Emax_C3aRblock*power(C_pep,gamma_pep)/(power(EC50_pep,gamma_pep)+power(C_pep,gamma_pep))'
 'Ps'          , 'Emax_switch   *power(C_pep,gamma_pep)/(power(EC50_pep,gamma_pep)+power(C_pep,gamma_pep))'
 'Pm'          , 'Emax_migration*power(C_pep,gamma_pep)/(power(EC50_pep,gamma_pep)+power(C_pep,gamma_pep))'
 'prot'        , 'NTF/(EC50_ntf+NTF)'
 'rpe_d'       , 'k_rpe_stress*Stress + k_rpe_cyt*Cyt_pro + k_rpe_phago*M_mig'
 'rgc_d'       , '(k_rgc_cyt*Cyt_pro + k_rgc_rpe*(1-RPE) + k_rgc_iop*Stress)*(1-prot)'
 'mig'         , '(k_mig_damp*DAMPs + k_mig_C5a*C5a)*(1-Pm)*M0'
 'inhib_anti'  , '1/(1 + Cyt_anti/K_anti_M1)'
 'tlr4'        , 'k_damp_M1*DAMPs*M_mig*inhib_anti'
 'c3ar'        , 'k_C3aR_act*C3a*(1-R_des)*M_mig*(1-Pc)*inhib_anti'
 'M1sw'        , 'k_M1_switch*(1+Ps)*M1'
 'C3cl'        , 'k_C3_cleave*C1q*C3' };
for i = 1:size(inter,1)
    p = addparameter(m, inter{i,1}, 0);
    p.ConstantValue = false;                       % required for assignment target
    addrule(m, sprintf('%s = %s', inter{i,1}, inter{i,2}), 'repeatedAssignment');
end

% ========================================================================
% STEP 5.  The 17 rate rules  (dX/dt = RHS), term-for-term from the R model.
% ========================================================================
rate = {
 'RPE'      , '-rpe_d*RPE'
 'DAMPs'    , 'k_damp_rpe*rpe_d*RPE + k_damp_rgc*rgc_d*RGC - k_damp_clear*DAMPs'
 'M0'       , '-mig + k_deact_M1*M1 + k_res_M2*M2 + k_return*M_mig'
 'M_mig'    , 'mig - c3ar - tlr4 - k_return*M_mig'
 'M1'       , 'tlr4 + c3ar - M1sw - k_deact_M1*M1'
 'M2'       , 'M1sw - k_res_M2*M2'
 'C1q'      , 'k_C1q_M1*M1*(1-Pb) - k_C1q_deg*C1q'
 'C3'       , 'k_C3_base + k_C3_rpe*rpe_d*RPE - C3cl - k_C3_deg*C3'
 'C3a'      , 'k_C3a_frac*C3cl - k_C3a_deg*C3a'
 'C5a'      , 'k_C5a_frac*C3cl - k_C5a_deg*C5a'
 'Cyt_pro'  , 'k_M1_cyt*M1 - k_deg_pro*Cyt_pro - k_inhib*Cyt_anti*Cyt_pro'
 'Cyt_anti' , 'k_M2_cyt*M2 - k_deg_anti*Cyt_anti'
 'NTF'      , 'k_ntf_base + k_M2_ntf*M2 - k_deg_ntf*NTF'
 'RGC'      , '-rgc_d*RGC'
 'A_eye'    , '-k_abs*A_eye'
 'C_pep'    , 'k_abs*A_eye - k_el_pep*C_pep'
 'R_des'    , 'k_des_on*C3a*(1-R_des) - k_des_off*R_des' };
for i = 1:size(rate,1)
    addrule(m, sprintf('%s = %s', rate{i,1}, rate{i,2}), 'rate');
end

% ========================================================================
% STEP 6.  Dosing.  *** CHOOSE THE SCHEDULE DELIBERATELY ***
%   doseTimes = [0 90 150 210];  % dashboard-faithful (4 admins)
%   doseTimes = [0 150 210];     % R-file-faithful (day-90 dropped)
%   doseTimes = [90 150 210];    % literal dose_times intent (no loading dose)
% ========================================================================
doseTimes  = [0 90 150 210];
doseAmount = 5.0;
dose = sbiodose('pep_dose', 'schedule');
dose.TargetName = 'A_eye';
dose.Time       = doseTimes;
dose.Amount     = repmat(doseAmount, size(doseTimes));

% ========================================================================
% STEP 7.  Solver — ode15s (stiff), tight tolerances to mirror lsoda.
% ========================================================================
cs.SolverType                        = 'ode15s';
cs.StopTime                          = 365;
cs.SolverOptions.AbsoluteTolerance   = 1e-10;
cs.SolverOptions.RelativeTolerance   = 1e-8;
cs.SolverOptions.MaxStep             = 1;        % cap step so dose times resolve
cs.SolverOptions.OutputTimes         = 0:1:365;  % daily output

% ========================================================================
% STEP 8.  Simulate treated + control, validate, plot.
% ========================================================================
if nargout == 0
    sdT = sbiosimulate(m, cs, [], dose);   % treated
    sdC = sbiosimulate(m, cs, [], []);     % control (no dose)

    get1 = @(sd,name) sd.Data(:, strcmp(sd.DataNames, name));
    t    = sdT.Time;
    fprintf('Treated RGC endpoint: %.3f%%\n', 100*get1(sdT,'RGC')(end));
    fprintf('Control RGC endpoint: %.3f%%\n', 100*get1(sdC,'RGC')(end));

    % VALIDATION: microglia pool invariant M0+M_mig+M1+M2 == M_total
    pool = get1(sdT,'M0')+get1(sdT,'M_mig')+get1(sdT,'M1')+get1(sdT,'M2');
    fprintf('Max |microglia pool - M_total| (treated): %.2e\n', ...
            max(abs(pool - pv.M_total)));

    figure;
    subplot(2,1,1); hold on;
    plot(t, 100*get1(sdT,'RGC'), 'LineWidth',1.5);
    plot(t, 100*get1(sdC,'RGC'), '--', 'LineWidth',1.5);
    ylabel('RGC survival (%)'); legend('Treated','Control'); title('RGC');
    subplot(2,1,2); hold on;
    plot(t, get1(sdT,'M1'),  'LineWidth',1.5);
    plot(t, get1(sdT,'C3a'), 'LineWidth',1.5);
    plot(t, get1(sdT,'C1q'), 'LineWidth',1.5);
    ylabel('amount'); xlabel('day'); legend('M1','C3a','C1q');
    title('Complement–microglia loop (treated)');
end
end

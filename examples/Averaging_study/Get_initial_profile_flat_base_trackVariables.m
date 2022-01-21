%% Get a time-steady profile (ice-sheet and subglacial) from a good guess.
% First written 29/3/21
% Driver script to run full ice-sheet--subglacial hydrology model
% for some predefined amount of time.

clear all;
addpath('../src/');
figpref(4);

%-------------------------------------------------------------------------%
%               Much of this is redundant because I load the correct
%               parameters
%-------------------------------------------------------------------------%


%Define dimensionless parameters and solver specifications
%specs  = solverSpecs;
%
%specs.npoints = 159;
%specs.HSolver = 'finiteVolumeUpwind';
%specs.usolver = 'fsolveCentredDiff';
%
%%Choose optional settings 
%specs.calvingType     = 'topographyDependent';
%specs.helConstitutive = 'exponential';
%specs.msType          = 'msSweep';
%
%%Set parameter values In this case it doesn't matter because the parameters
%%are those of the initial guess.
%specs.realisticParams = 0;
%if specs.realisticParams ==1
%    %calculate dimensionless parameters based on realistic dimensional paramters
%    params = assignNonDim(dimparams);
%else
%    % Use a chosen (easy to work with) set of dimenionless parameters
%    %params = nonDimParams;
%    params = assignNonDim(dimparams);
%    
%    params.epsilon    = 0.001;
%    delta             = 0.001;
%    params.delta      = delta;
%    params.sigmatilde = delta*0.1;
%    params.n          = 3;
%    params.gamma      = 1;
%    params.lambdatilde= 0.1;
%    params.omega_r    = 10;
%    params.Atilde     = params.omega_r/params.delta/params.n^(params.n);
%    params.mutilde    = 1;
%    params.Ktilde     = 0.3;
%end


% Load initial guess
%guess_root = './Initial_guess/noHel/';
%initial_sol = load([guess_root 'LargestSigma_long_r2_rs6p66_n3.csv']);
%specs.sol_0 = initial_sol;
fn       = ['./data/flat_base/16_4_21_noHel/JustBelowSL_higher_melt_lambda0p1_eps0p001_n3deltams1.csv'];
initial_sol = load(fn);
% Load parameters corresponding to initial guess
params_root = './data/flat_base/16_4_21_noHel/';
pobj = load([params_root 'params_n3.mat']);
params = pobj.params;
sobj = load([params_root 'specs_n3.mat']);
specs = sobj.specs;

specs.sol_0 = initial_sol;

%-------------------------------------------------------------------------%
%                       Integrate in time
%-------------------------------------------------------------------------%
params.period = 10;
specs.t_max   = params.period; %Note that this may be insufficient if the initial guess isn't good enough
specs.tpoints = 5000;
specs.increments = 50;
specs.advanceGlacier = 1;
specs.advanceSubglacial = 1;
specs.msAmplitude  = 0; %No melt rate osscilation 

% Change a parameter from initial guess
%params.rs = params.r/0.39;
params.b = params.b_xvals*0 - 0.05;

params.sigmatilde = 1;
%params.rs = params.r/0.22;


%Set meltrate based on initial profile
params.ms          = surfaceMelt(specs.sol_0, params, specs, 0, params.ms);
params.ms_subOnly  = surfaceMelt(specs.sol_0, params, specs, 0, params.ms);


% Run the solver until stationary or until t_{max} is reached.
[sol, params, specs, trackedVariables]  = run_solver_Lfixed_trackVariables(params, specs, 1);

writematrix(trackedVariables.Hcalving, './data/Hcalving_long.csv');
writematrix(trackedVariables.dSol, './data/dSol_long.csv');

out_root = './Initial_guess/noHel/';
fn       = [out_root 'LargeSigma_long_r2_rs6p66_n' int2str(params.n) '.csv'];
writematrix(sol, fn)

save([out_root 'paramsLong_n3.mat'], 'params')
save([out_root 'specsLong_n3.mat'], 'specs')


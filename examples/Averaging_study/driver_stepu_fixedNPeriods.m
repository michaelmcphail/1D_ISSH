%% Driver for running subglacial--ice-sheet-momentum equations for a fixed number of periods
clear all;
addpath('../../src/');
figpref(4);

%Define dimensionless parameters and solver specifications
   
%                           Load initial state
%-------------------------------------------------------------------------%

%Note: be careful about load order, only set initial state after loading
%specs because other wise you will simply overwrite prescribed initial state

fn       = ['./data/flat_base/16_4_21_noHel/JustBelowSL_higher_melt_lambda0p1_eps0p001_n3deltams1.csv'];

params_root = './data/flat_base/16_4_21_noHel/';
pobj = load([params_root 'params_n3.mat']);
params = pobj.params;
sobj = load([params_root 'specs_n3.mat']);
specs = sobj.specs;

%Load initial state
initial_steadysol =load(fn);
specs.sol_0 = initial_steadysol(:, 1);


%-------------------------------------------------------------------------%
%                Integrate in time (sweep through ms amplitudes)
%-------------------------------------------------------------------------%
params.period = 0.001;
specs.t_max   = params.period;
specs.tpoints = 70000;
specs.increments = 200;
specs.advanceGlacier = 0;
specs.advanceSubglacial = 1;

deltamsList = [0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];

specs.periodicity_tol = 0.0005;

specs.msAmplitude  = deltamsList(1);



% Run fully coupled for post-averaging solution
nupdate_list = [0];



for k = 1:length(nupdate_list)
    n_updates = nupdate_list(k);
    root = ['./data/flat_base/3_11_21_noHel_stepu' int2str(n_updates) '/'];
    
    fnFull_root =[root 'stepu_nUpdates'];
    fntvals_root = [root 'tvals_stepu_nUpdates'];
    try
        % Use velocity calculated from average to advance N (i.e. the coarse
        % method)

          for i = 6
            specs.msAmplitude = deltamsList(i);
            [sol, params, specs] = find_periodic_fixed_sheet_stepu_fixedNPeriods(params, specs, n_updates);
            specs.sol_0 = initial_steadysol(:, 1);
    
            if ~exist(root, 'dir')
                mkdir(root)
            end
            fieldsfn       = [ fnFull_root int2str(n_updates) 'deltams' int2str(i) '.csv'];
            writematrix(sol, fieldsfn)
            
            
            tvals = linspace(0, specs.t_max, specs.increments + 1);
            tvalsfn = [fntvals_root int2str(n_updates) 'deltams' int2str(i) '.csv'];
            writematrix(tvals, tvalsfn)
            
        end
        
    catch ME
        fprintf('coarse solution without success: %s\n', ME.message);       
    end
    
end

save([root 'params_n3.mat'], 'params')

save([root 'specs_n3.mat'], 'specs')

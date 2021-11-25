function [sol_out, paramsOut, specsOut] = find_periodic_fixed_sheet_stepu_fixedNPeriods(params, specs, n_updates)
%  FIND_PERIODIC_FIXED_SHEET_STEPU_FIXEDNPERIODS run subglacial model,
%  updating the momentum equation n_updates times per period, until
%  difference between periods is below tolerance.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.
%    n_updates: number of velocity updates to perform within a period.
%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%
    npointsOld = specs.npoints;                            % Number of z points
    t          = 0;                                        % Initial time
    dtmax      = specs.t_max/specs.tpoints;                % Time step
    z_old      = linspace(0, specs.z_max, specs.npoints);  % grid for old solution

    epsilon   = params.epsilon;

    % first save values to save output vector
    Hold     = specs.sol_0(1:npointsOld);
    uold     = specs.sol_0(npointsOld+1:2*npointsOld);
    Nold     = specs.sol_0(2*npointsOld+1:3*npointsOld);
    phiold   = specs.sol_0(3*npointsOld+1:4*npointsOld);
    h_old    = specs.sol_0(4*npointsOld+1:5*npointsOld);
    helold   = specs.sol_0(5*npointsOld+1:6*npointsOld);
    hcavold  = specs.sol_0(6*npointsOld+1:7*npointsOld);
    Lold     = [specs.sol_0(end)];
    
    
    %Specify when to save outputs (checkpoints)
    chkPointdeltat = specs.t_max/specs.increments;
    sol_out        = specs.sol_0;
    chkPointNum    = 1;
    




    %New specs based on new epsilon
    %epsilon
    npoints    = min(max(npointsOld, ceil(5/sqrt(epsilon))), 300); %resolution based on velocity boundary layer (5 points in boundary layer))
    dz         = 1/npoints;
    dt         = min([0.9*dz/max(uold), dz/5, 0.001]); % Time step satisfying an estimated CFL condition;


    disp(['For epsilon = '  num2str(epsilon) ', we are using ' num2str(npoints) ' gridpoints and the time step is ' num2str(dt) '.'])

    specs.npoints = npoints;
    z_values      = linspace(0, specs.z_max, specs.npoints);

    %Interpolate old solution onto new grid as the initial guess
    Hinterp     = interp1(z_old', Hold, z_values');
    uinterp     = interp1(z_old', uold, z_values');
    Ninterp     = interp1(z_old', Nold, z_values');
    phiinterp   = interp1(z_old', phiold, z_values');
    hinterp     = interp1(z_old', h_old, z_values');
    helinterp   = interp1(z_old', helold, z_values');
    hcavinterp  = interp1(z_old', hcavold, z_values');



    % Initial state
    sol = [Hinterp; uinterp; Ninterp; phiinterp; hinterp; helinterp; hcavinterp; Lold];

    Hvals     = Hinterp;
    uvals     = uinterp;
    Nvals     = Ninterp;
    phivals   = phiinterp;
    hvals     = hinterp;
    helvals   = helinterp;
    h_cavvals = hcavinterp;
    Lval      = Lold;

%-------------------------------------------------------------------------%
%                   Integrate in time
%-------------------------------------------------------------------------%
    count      = 1;
    msOriginal = params.ms;
    nInternalIts = 1;
    dt = min(dtmax, 0.5*dz/max(abs(uvals)));
    
    ntpoints          = specs.increments; %check that this is the number of tpts
    solOld            = sol*ones(1, ntpoints); %check dimensions
    periodicity_error = 1;
    NPeriods = 1;
    for k = 1:20
        % Advance solution through one period
        [sol] = run_subglacial_one_period_stepu_ctsAveraging(params, specs, n_updates);
        %check how periodic solution is (by comparing to previous period)
        if length(solOld(1, :)) == length(sol(1, :))
            periodicity_error = norm(solOld - sol)/ntpoints/length(sol(:, 1));
        else
            periodicity_error = 1;
        end
        specs.sol_0 = sol(:, end);
        
        solOld = sol;
        
        fprintf(' periodicity error = %f, number of periods = %i \n', periodicity_error, count);
        count = count + 1;
        
        %Save Diagnostic variables to temp the last successful period
        if specs.msAmplitude >0.8
            diagnosticRoot = './data/temp/3_11_21_stepu_20periodLarge/';
        else
            diagnosticRoot = './data/temp/3_11_21_stepu_20periodSmall/';
        end
        if ~exist(diagnosticRoot, 'dir')
            mkdir(diagnosticRoot)
        end
        fieldsfn       = [ diagnosticRoot 'updates' int2str(n_updates) '.csv'];
        
        writematrix(sol, fieldsfn);
        fid = fopen([diagnosticRoot 'errOut.txt'], 'a+');
        fprintf(fid, '%i %f %f\n', [n_updates, specs.msAmplitude, periodicity_error]);
        fclose(fid);
        
        tvals = linspace(0, specs.t_max, specs.increments + 1);
        tvalsfn = [diagnosticRoot 'tvals_updates' int2str(n_updates) '.csv'];
        writematrix(tvals, tvalsfn)



    end
    sol_out  = sol;
    
    

    paramsOut = params;
    specsOut  = specs;

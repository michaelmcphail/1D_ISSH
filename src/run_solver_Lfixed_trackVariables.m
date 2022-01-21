function [sol_out, paramsOut, specsOut, trackedVariables] = run_solver_Lfixed_trackVariables(params, specs, sm_continuation_param)
%  RUN_SOLVER_LFIXED_TRACKVARIABLES run fully coupled model (ice sheet+
%  subglacial) with the end point pinned, for some fixed interval of time.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.
%    sm_continuation_param: scale factor of surface melt, useful for starting
%        from zero surface melt.

%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%
    npointsOld = specs.npoints;                            % Number of z points
    t          = 0;                                        % Initial time
    %dtmax      = specs.t_max/specs.tpoints;                % Time step
    dtmax      = 0.0001;
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

    trackedVariables = struct;

    trackedVariables.Hcalving = [];
    trackedVariables.dSol = [];

%-------------------------------------------------------------------------%
%                   Integrate in time
%-------------------------------------------------------------------------%
    t_max       = specs.t_max;
    %t_max       = 0.5;
    count      = 1;
    msOriginal = params.ms;
    dSol = 1;
    while t<t_max && dSol > 0.0001
        dt = min(dtmax, 0.5*dz/max(abs(uvals)));
        solOld = sol;
        % Time step the thickness $H$ by dt
        sol  = timestep_Thickness_fixedL(sol, dt, params, specs);


        % Prescribe the number of sub-iterations for the coupled momentum
        % equation--subglacial system (minimum of 30, arbitrary).
        nInternalIts = max(ceil(dt/min([0.5*params.delta])) , 30);

        for i = 1:nInternalIts
            %Increment by dt/nInternalIts
            % Should do u subscycle?
            sol  = solve_Velocity(sol, params, specs);
            sol  = timestep_Subglacial_helConstitutive_fixedL(sol, dt/nInternalIts, params, specs);
        end

        if t > chkPointNum*chkPointdeltat
            %interpolate between sol(t+dt) and old_sol(t), add to sol_out
            %For now just use weighted averate8
            tchk = chkPointNum*chkPointdeltat;
            t_old = t - dt;
            sol_incr = (solOld*(tchk - t_old) + sol*(t - tchk))/dt;
            chkPointNum = chkPointNum + 1;
        elseif t == chkPointNum*chkPointdeltat
           chkPointNum = chkPointNum + 1;

        end
        
        Hvals     = sol(1:npoints);
        uvals     = sol(npoints+1:2*npoints);
        Nvals     = sol(2*npoints+1:3*npoints);
        phivals   = sol(3*npoints+1:4*npoints);
        hvals     = sol(4*npoints+1:5*npoints);
        helvals   = sol(5*npoints+1:6*npoints);
        h_cavvals = sol(6*npoints+1:7*npoints);
        Lval      = sol(end);


        dSol = norm((solOld -sol)/dt/npoints);
        fprintf('t = %f, dSol = %f, dt = %f, Hval = %f, count = %i \n', t, dSol, dt, Hvals(end), count);

        t =t + dt;
        count = count + 1;

        %Get new surface melt rate
        msOut             = sm_continuation_param*surfaceMelt(sol, params, specs, 0, msOriginal);
        params.ms         = msOut;
        params.ms_subOnly = msOut;
        
        trackedVariables.Hcalving = [trackedVariables.Hcalving, Hvals(end)];
        trackedVariables.dSol = [trackedVariables.dSol, dSol];


    end

    %store final solution
    sol_out  = sol;
    paramsOut = params;
    specsOut  = specs;

function [sol_out] = run_subglacial_one_period_stepu_ctsAveraging(params, specs, n_updates)
%  RUN_SUBGLACIAL_ONE_PERIOD_STEPY_CTSAVERAGING run the subglacial model,
%  updating the average effective pressure continuously, for one period
%  with n_updates throughout.
%    params: contains the dimensionless parameters
%    specs: contains the solver specifications
%    n_updates: number of velocity updates to perform within a period

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
    t_max       = params.period;
    msOriginal  = params.ms;
    nInternalIts = 1;    
    %Specify when to save outputs (checkpoints)
    chkPointdeltat = t_max/specs.increments;
    chkPointNum    = 1;
    count = 1;
    tvals_checkpoints = linspace(0, t_max, specs.increments);
    
    updateInterval = t_max/(n_updates+1);
    updateNumber   = 1;
    
    sol_ave = sol;
    time_since_update = 0;

    while t<t_max
        solOld = sol;
        dt = min(dtmax, 0.5*dz/max(abs(uvals)));
        
        sol  = timestep_Subglacial_helConstitutive(sol, dt, params, specs);

        sol_ave = (time_since_update*sol_ave + dt*sol)/(time_since_update + dt);
        time_since_update = time_since_update + dt;
        
        %Do checkpoint saving here!
        if t > chkPointNum*chkPointdeltat
            %interpolate between sol(t+dt) and old_sol(t), add to sol_out
            %For now just use weighted averate8
            tchk = chkPointNum*chkPointdeltat;
            t_old = t - dt;
            
            sol_incr = (solOld*(tchk - t_old) + sol*(t - tchk))/dt;
            sol_out  = [sol_out, sol_incr];
            
            
            chkPointNum = chkPointNum + 1;
        elseif t == chkPointNum*chkPointdeltat
            sol_out  = [sol_out, sol];
            chkPointNum = chkPointNum + 1;
            
        end

        %If on or passed update time then update the velocity
        if t > updateNumber*updateInterval
            %Insert updated velocity into vector
 
            %TODO: interpolate to account for discrepancy between current
            %time and update time (as per below)
            %interpolate to get solution at tupdate
            %tupdate = updateNumber*updateInterval;
            %t_old = t - dt;
            %sol_incr = (solOld*(tupdate - t_old) + sol*(t - tupdate))/dt;            
            %sol_ave = ((time_since_update - dt)*sol_ave + sol_incr*(tupdate - (time_since_update - dt)))/updateInterval;

            %Reset time since velocity update
            time_since_update = 0;
            
            try
            sol_updatedVel  = solve_Velocity(sol_ave, params, specs);
            catch MESSAGE
                fakeSet = 1;
               disp('oh no!'); 
            end

            %Insert predicted velocity
            sol = [sol(1:specs.npoints); sol_updatedVel((specs.npoints+1):2*specs.npoints); sol((2*specs.npoints+1):end)];
            updateNumber = updateNumber + 1;
       
        elseif t == updateNumber*updateInterval
            %Get averaged N --- Not quite right, would probably need to account for residuals at either end, will probably end up using even intervals; in which case it doesn't matter.
            %tupdate = updateNumber*updateInterval;


            %sol_ave =sol_out(tvals_checkpoints>(updateNumber-1)*updateInterval&tvals_checkpoints <tupdate);
            %sol_ave = sum(sol_ave, 2)/length(sol_ave(1, :));
            sol_updatedVel  = solve_Velocity(sol_ave, params, specs);


            sol = [sol(1:specs.npoints); sol_updatedVel((specs.npoints+1):2*specs.npoints); sol((2*specs.npoints+1):end)];
            time_since_update = 0;
            updateNumber = updateNumber + 1;
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
        fprintf('t = %f, dSol = %f, dt = %f, count = %i, melt = %f \n', t, dSol, dt, count, norm(surfaceMelt(sol, params, specs, t, msOriginal)));
        
        t =t + dt;
        count = count + 1;
        
        %Get new surface melt rate
        msOut             = surfaceMelt(sol, params, specs, t, msOriginal);
        params.ms         = msOut;
        params.ms_subOnly = msOut;




    end
    %If while loop finishes before updating, make last update.
    if time_since_update>0
        sol_updateVel  = solve_Velocity(sol_ave, params, specs);
        sol = [sol(1:specs.npoints); sol_updatedVel((specs.npoints+1):2*specs.npoints); sol((2*specs.npoints+1):end)];
    end
    %store final solution
    sol_out  = [sol_out, sol];

    paramsOut = params;
    specsOut  = specs;

function [sol] = timestep_Thickness_fixedL(sol, dt, params, specs)
%  TIMESTEP_THICKNESS_FIXEDL step the ice-sheet thickness forwards in time
%  by dt
%    sol: full solution vector at current time.
%    dt: time step to advance solution by.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.
%    sm_continuation_param: scale factor of surface melt, useful for starting
%        from zero surface melt.

%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%

    %Specifications
    npoints     = specs.npoints;
    dz          = specs.z_max/npoints;
    z_values    = linspace(0, specs.z_max, npoints);
    t_values    = linspace(0, specs.t_max, specs.tpoints);
    regparam    = specs.regparam; % regularisation term for viscosity

    %Prescribed functions and exponents
     %Defined profiles.
     b_xvals     = params.b_xvals;    % x coordinates of basal data
     a_xvals     = params.a_xvals;    % x coordinates of accumulation data
     ms_xvals    = params.ms_xvals;   % x coordinates of surface melt data
     bvals_atx   = params.b;          % Basal topography
     avals_atx   = params.a;          % Accumulation rate
     msvals_atx  = params.ms;         % Surface melt rate
     mbtilde     = params.mbtilde;    % Basal melt rate

    alpha   = params.alpha;  % Hydraulic transmissivity exponent
    n       = params.n;      % Glen's law coefficient
    gamma   = params.gamma;  % h_el exponent

    %Relevant dimensionless parameters
     %Ice properties (momentum equation)
    rhotilde_w  = params.rhotilde_w;  % Ratio of water to ice density.
    mutilde     = params.mutilde;     % Dimensionless basal friction coeff.
    beta        = params.beta;        % Dimensionless side friction coeff.
    lambdatilde = params.lambdatilde; % Dimensionless ~roughness.
    epsilon     = params.epsilon;     % Dimensionless extensional viscosity coeff.
    rhotilde_o  = params.rhotilde_o;  % Ratio of ocean to fresh water density.

     %Subglacial parameters
    Ktilde     = params.Ktilde;       % Dimensionless hydraulic conductivity
    sigmatilde = params.sigmatilde;   % Dimensionless rate of englacial storage change
    delta      = params.delta;        % Dimensionless rate of sheet height change
    nu         = params.nu;           % Ratio of Nscale to phiscale
    omega_r    = params.omega_r;      % ~Dimensionless rate of cavitation
    Atilde     = params.Atilde;       % Dimensionless rate of viscous
    hrtilde    = params.hrtilde;      % Dimensionless basal roughness
    hctilde    = params.hctilde;      % Dimensionless coefficient of h_el


    %Fornberg weights for order dz^3 accuracy (curently unused)
    fornberg_weights = [-11/6, 3, -3/2, 1/3];

    % value of solution at $t$
    Hvals     = sol(1:npoints);
    uvals     = sol(npoints+1:2*npoints);
    Nvals     = sol(2*npoints+1:3*npoints);
    phivals   = sol(3*npoints+1:4*npoints);
    hvals     = sol(4*npoints+1:5*npoints);
    helvals   = sol(5*npoints+1:6*npoints);
    h_cavvals = sol(6*npoints+1:7*npoints);
    Lval      = sol(end);

    %Derived values
    bvals_z   = interp1(b_xvals, bvals_atx, Lval*z_values');   % Basal topography (at z points)
    avals_z   = interp1(a_xvals, avals_atx, Lval*z_values');   % Accumulation rate (at z points)
    msvals_z  = interp1(ms_xvals, msvals_atx, Lval*z_values'); % Surface melt rate (at z points)

    %-------------------------------------------------------------------------%
    %                   Update H
    %-------------------------------------------------------------------------%


    switch specs.HSolver
      case 'forwardsEuler'
       % Update $H$ using the forwards Euler method (in conjunction with method of lines).
          dH0   = - Hvals(1).*(uvals(2) - uvals(1))/dz/Lval + avals_z(1) - msvals_z(1) - mbtilde;
          dH    = -(uvals(2:end)).*(Hvals(2:end) - Hvals(1:end-1))/dz/Lval - Hvals(2:end).*(uvals(2:end) - uvals(1:end-1))/dz/Lval +avals_z(2:end) - msvals_z(2:end) - mbtilde;

          Hvals = Hvals + dt*[dH0; dH];
          Lval  = Lval ;
      case 'backwardsEuler'
       % Update $H$ using the backwards Euler method (in conjunction with method of lines).
          [newsol, fval,exitflag,output] = fsolve(@backwardsEuler, [Hvals; Lval]);

          if exitflag<0
              disp('Error: Solution not found')
              return;
          end

          Hvals = newsol(1:end-1);
          Lval  = newsol(end);


       case 'finiteVolumeUpwind'
       % A finite volume method emulating BISICLES (Godunov method. Because of flux integral approximation this is now equivalent to UPWIND I beleive)
           Hmidpoints    =  Hvals(1:end-1); % Assuming velocity is positive
           umidpoints    = (uvals(2:end) + uvals(1:end-1))/2;
           fluxmidpoints = Hmidpoints.*umidpoints/Lval ;
           sourceterm    = avals_z - msvals_z - mbtilde;
           fluxend       = Hvals(end)*uvals(end)/Lval;


           dH0   = -2*fluxmidpoints(1)/dz +sourceterm(1);
           dH    = (fluxmidpoints(1:end-1) - fluxmidpoints(2:end))/dz + sourceterm(2:end-1);
           dHn   = (fluxmidpoints(end) - fluxend)/dz + sourceterm(end); %TODO figure out what to actually do with ghost point here!

           Hvals = Hvals + dt*[dH0; dH; dHn];
           Lval  = Lval;

        case 'finiteVolumeFaceCentred'
        % A finite volume method using a cell-centred--face-centred approach
            Hmidpoints    = (Hvals(2:end) + Hvals(1:end-1))/2;
            umidpoints    = (uvals(2:end) + uvals(1:end-1))/2;
            fluxmidpoints = Hmidpoints.*umidpoints/Lval ;
            sourceterm    = avals_z - msvals_z - mbtilde ;
            fluxend       = Hvals(end)*uvals(end)/Lval;


            dH0   = -2*fluxmidpoints(1)/dz +sourceterm(1);
            dH    = (fluxmidpoints(1:end-1) - fluxmidpoints(2:end))/dz + sourceterm(2:end-1);
            dHn   = (fluxmidpoints(end) - fluxend)/dz + sourceterm(end); %TODO figure out what to actually do with ghost point here!

            Hvals = Hvals + dt*[dH0; dH; dHn];
            Lval  = Lval;
    end

    %Return sol with updated Hvals
    sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];




    %-------------------------------------------------------------------------%
    %                   Backwards Euler objective function
    %-------------------------------------------------------------------------%

        function res = backwardsEuler(solve0)
            Hnew = solve0(1:end-1);
            Lnew = solve0(end);

            %For fully implicit, calvingRate should take in appropriate solution, not sure how much of a difference it will make.
            dH0 = - Hnew(1).*(uvals(2) - uvals(1))/dz/Lval + avals_z(1) - msvals_z(1) - mbtilde;
            dH  = -(uvals(2:end) ).*(Hnew(2:end) - Hnew(1:end-1))/dz/Lnew - Hnew(2:end).*(uvals(2:end) - uvals(1:end-1))/dz/Lnew +avals_z(2:end) - msvals_z(2:end) - mbtilde;


            resdL  = Lnew  - Lval;
            resdH0 = Hnew(1) - Hvals(1) - dt*dH0;
            resdH  = Hnew(2:end) - Hvals(2:end) - dt*dH;


            %Residual vector
            res    = [resdH0; resdH; resdL];

        end

    %-------------------------------------------------------------------------%
    %                   Other functions
    %-------------------------------------------------------------------------%
    % Event function to update initial guesses (if using Matlab integrator)
        % function f = events(z, n)
        %     H_at_z = n(1:end-1);
        %     L_at_z = n(end);
        %     phi_at_t = shoot_for_phi(specs.phiguess, Hvals, Lval, params, specs); %TODO Write this
        %     specs.phiguess = phi_at_t(1);
        %     Nvals_at_t   = gravPotVals - phivals;
        %     [~, newv0] = shoot_for_u_oddn(specs.v0guess, Nvals_at_t, H_at_z, L_at_z, params, specs);
        %     specs.v0guess = newv0;
        %
        % end


    end

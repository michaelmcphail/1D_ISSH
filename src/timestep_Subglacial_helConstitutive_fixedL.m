function [sol] = timestep_Subglacial_helConstitutive_fixedL(sol, dt, params, specs)
%  TIMESTEP_SUBGLACIAL_HELCONSTITUTIVE advance the subglacial system forward 
%  in time by dt
%    sol: full solution vector at current time.
%    dt: time step to advance solution by.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.
%    sm_continuation_param: scale factor of surface melt, useful for starting


%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%

    %Specifications
    npoints     = specs.npoints;
    dz          = specs.z_max/npoints;
    z_values    = linspace(0, specs.z_max, npoints);
    t_values    = linspace(0, specs.t_max, specs.tpoints);
    regparam    = specs.regparam; % regularisation term for viscosity
    JFNKtol     = specs.JFNKtol; %Tolerance for Newton iterations

    %Prescribed functions and exponents
     %Defined profiles.
    b_xvals     = params.b_xvals;    % x coordinates of basal data
    a_xvals     = params.a_xvals;    % x coordinates of accumulation data
    ms_xvals    = params.ms_xvals;   % x coordinates of surface melt data
    bvals_atx   = params.b;          % Basal topography
    avals_atx   = params.a;          % Accumulation rate
    msvals_atx  = params.ms_subOnly;         % Surface melt rate
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


    % Derived values
    bvals_z   = interp1(b_xvals, bvals_atx, Lval*z_values');   % Basal topography (at z points)
    avals_z   = interp1(a_xvals, avals_atx, Lval*z_values');   % Accumulation rate (at z points)
    msvals_z  = interp1(ms_xvals, msvals_atx, Lval*z_values'); % Surface melt rate (at z points)

    midpoints   = z_values(1:end-1) + z_values(2)/2;    % z_values in between grid points
    b_end       = bvals_z(end);               % topography height at b_end


%-------------------------------------------------------------------------%
%                   Update subglacial hydrology
%-------------------------------------------------------------------------%

    switch specs.subglacialSolver
        case 'centredDiff'
            phivals_old     = phivals;
            h_cavvals_old = h_cavvals;
            options = optimset('Display', 'off', 'TolFun', 10^(-12));
            [w,fval,exitflag,output]  = fsolve(@finite_diff_phi_hcav, [phivals; h_cavvals], options);

            if exitflag<0
                disp('Solution not found')
                return;
            end

            phivals     = w(1:npoints);
            h_cavvals = w(npoints+1: 2*npoints);
            Nvals   = (Hvals + rhotilde_w*bvals_z - phivals)/nu;
            [helvals, ~] = hel_Constitutive(sol, params, specs);
            hvals     = helvals + h_cavvals;

        case 'JFNK'
            %Jacobian-Free Newton Krylov method using Picard linearisation on differential operator
             nResid  = JFNKtol + 1;
             JFNKits = 0; %current number of JFNK iterations
             while nResid > JFNKtol && JFNKits<specs.JFNKitsMax
                   [Jac, Rhs]   = linearised_subglacial_op(sol, dt, params, specs);
                   linearRes    = Jac*phivals - Rhs;
                   dphi         = Jac\(-linearRes);
                   [w, res, fl] = linesearch_w_phi(sol, dt, params, specs, dphi, linearRes);
                   if fl == 1
                     phivals   = phivals + w*dphi;
                     Nvals     = (Hvals + rhotilde_w*bvals_z - phivals)/nu;
                     h_cavvals = update_h_cav(h_cavvals, Nvals, dt);
                     [helvals, ~] = hel_Constitutive(sol, params, specs);
                     hvals     = helvals + h_cavvals;
                   else
                     phivals = phivals;
                     disp('Oh No!')
                     break;
                   end
                   sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
                   nResid   = norm(res);
                   JFNKits  =JFNKits+1;
            end
      end

      sol = [Hvals; uvals;Nvals;phivals; hvals; helvals; h_cavvals; Lval];


%-------------------------------------------------------------------------%
%                   Residual functions
%-------------------------------------------------------------------------%

% Centred difference for mass conservation and backwards Euler for h_cav evolution
    function res = finite_diff_phi_hcav(w)
        %Find an updated N and h_cav for a given H and u
        phivals     = w(1:npoints);
        h_cavvals = w(npoints+1: 2*npoints);
        
        Nvals   = (Hvals + rhotilde_w*bvals_z - phivals)/nu;
        [helvals, fprimeN] = hel_Constitutive(sol, params, specs);
        hvals     = helvals + h_cavvals;
        
        h_midpoints = interp1(z_values, hvals, midpoints)';
        
        %z-derivatives of hcav and N
        dhcavdz = [(h_cavvals(2) - h_cavvals(1))/dz; (h_cavvals(2:end) - h_cavvals(1:end-1))/dz];
        dphidz    = [(phivals(2) - phivals(1))/dz; (phivals(2:end) - phivals(1:end-1))/dz];
        
        
        dHdt    = avals_z(2:end) - msvals_z(2:end) - mbtilde - (Hvals(2:end).*uvals(2:end) - Hvals(1:end-1).*uvals(1:end-1))/dz;
        artificialDiffusion = specs.artificialDiffusitivity*([h_cavvals(2:end); h_cavvals(end-1)]-2*h_cavvals(1:end)+[h_cavvals(2); h_cavvals(1:end-1)] )/dz^2/Lval^2;
        dhcavdt = (rhotilde_w*mbtilde + omega_r*uvals.*max((1 - h_cavvals/hrtilde), 0) -delta* Atilde*sqrt(Nvals.^2 + regparam^2).^(n-1).*(max(Nvals, 0)).*h_cavvals )/delta + artificialDiffusion ;
        
        
        meltsource    = msvals_z(2:end -1) + mbtilde;
        nonlinLaplace = Ktilde*(h_midpoints(2:end).^3 .*(phivals(3:end) -phivals(2:end-1)) - h_midpoints(1:end-1).^3 .*(phivals(2:end-1) -phivals(1:end-2)))/dz^2/Lval^2;
        sheetchange   = delta*fprimeN(2:end -1)/nu.*dHdt(1:end-1);
        cavitation    = delta*dhcavdt(2:end-1);
        dphidtcoeff   = sigmatilde - delta*fprimeN(2:end-1)/nu;
        
        dphidt   = (meltsource +  nonlinLaplace - sheetchange - cavitation );
        
        
        
        
        
        
        resphi1      = phivals(2) - phivals(1);
        resphi2tonm1 = (phivals(2:end-1) - phivals_old(2:end-1)) - dt*dphidt./dphidtcoeff;
        %resNend    = -(rhotilde_w - rhotilde_o)*b_end*(b_end<0) + Hvals(end) + rhotilde_w*b_end - nu*Nvals(end);
        resphiend    = -(rhotilde_w - rhotilde_o)*b_end*(b_end<0) +phivals(end);
        reshcav    = h_cavvals - h_cavvals_old - dt*dhcavdt;
        
        
        
        res = [resphi1; resphi2tonm1; resphiend; reshcav];
        
        
    end

% Update h_cav given N using backwadsEuler

    function h_cav_updated = update_h_cav(h_cav_old, Nvals, dt)
        
        h_cav_updated = [(h_cav_old(1) + dt*rhotilde_w*mbtilde/delta)/(1 + dt*Atilde*sqrt(Nvals(1)^2 +regparam^2)^(n-1)*Nvals(1))];
        
        for i = 2:length(h_cav_old)
            a1 = rhotilde_w *mbtilde + omega_r*uvals(i);
            a2 = + omega_r*uvals(i)/hrtilde + Atilde*delta*sqrt(Nvals(i)^2 + regparam^2)^(n-1)*Nvals(i);
            new_h_cav_val = (h_cav_old(i) + dt*a1/delta)/(1 + dt*a2/delta);
            h_cav_updated = [h_cav_updated; new_h_cav_val];
        end
        
    end



end

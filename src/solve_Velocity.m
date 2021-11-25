function [sol] = solve_Velocity(sol, params, specs)
% SOLVE_VELOCITY solve momentum equation with $H$ and $N$ in sol. 
%    sol: full solution vector at current time.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.

%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%

    %Specifications
    npoints     = specs.npoints;
    dz          = specs.z_max/npoints;
    z_values    = linspace(0, specs.z_max, npoints);
    x_extended  = linspace(0, specs.x_extended_max, specs.nx_points);
    t_values    = linspace(0, specs.t_max, specs.tpoints);
    regparam    = specs.regparam; % regularisation term for viscosity
    JFNKtol     = specs.JFNKtol; %Tolerance for Newton iterations
    
    %Prescribed functions and exponents
    %Defined profiles; defined with respect to  x_extended.
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
    if isscalar(lambdatilde)
       lambdatilde = params.lambdatilde*sol(npoints+1:2*npoints).^0; 
    end
    epsilon     = params.epsilon;     % Dimensionless extensional viscosity coeff.
    rhotilde_o  = params.rhotilde_o;  % Ratio of ocean to fresh water density.

     %Subglacial parameters
    Ktilde     = params.Ktilde;       % Dimensionless hydraulic conductivity
    sigmatilde = params.sigmatilde;   % Dimensionless rate of change of englacial storage change
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
    H_midpoints = interp1(z_values, Hvals, midpoints)'; % Thickness in between gridpoints

    gravPotVals    = rhotilde_w*bvals_z + Hvals; % Gravitational potential values
    b_end          = bvals_z(end);               % topography height at b_end
    S_height       = bvals_z + Hvals;            % Surface height
    S_grad         = [(S_height(2) - S_height(1))/dz; (S_height(2:end) - S_height(1:end-1))/dz]; %not sure about initial grad

%-------------------------------------------------------------------------%
%                   Find velocity
%-------------------------------------------------------------------------%

    switch specs.usolver
        case 'fsolveCentredDiff'
          opts=  optimset('display','off');
            if n ==1
             % Newtonian fluid

                [uvals, fval,exitflag,output] = fsolve(@finite_diff_u_n1, uvals, opts);

                if exitflag<0
                    disp('Error: Solution not found')
                    return;
                end
            else
             % non-Newtonian fluid
                [uvals, fval,exitflag,output] = fsolve(@finite_diff_u, uvals, opts);

                if exitflag<0
                    disp('Error: Solution not found')
                    return;
                end
            end

        case 'JFNK'
            %Jacobian-Free Newton Krylov method using Picard linearisation on differential operator
             nResid  = JFNKtol + 1;
             JFNKits = 0; %current number of JFNK iterations
             while nResid > JFNKtol && JFNKits<specs.JFNKitsMax
                   [Jac, Rhs]   = linearised_mom_op(sol, params, specs);
                   linearRes    = Jac*uvals - Rhs;
                   du           = Jac\(-linearRes);
                   [w, res, fl] = linesearch_w(sol, params, specs, du, linearRes);
                   if fl == 1
                     %uvals = uvals + w*du;
                     uvals = uvals + du;
                   else
                     disp('could not reduce residual, using fsolve')
                       opts=  optimset('display','off');
                     %[uvals, fval,exitflag,output] = fsolve(@finite_diff_u_n1, uvals, opts);
                      uvals = uvals + du/32;
%                      if exitflag<0
%                          disp('Error: Solution not found')
%                          return;
%                      end
                     sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
                     %disp('Oh no! Could not find a suiable w using linesearch')
                     break;
                   end
                   sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
                   nResid   = norm(res)/npoints;
                   JFNKits  =JFNKits+1;
             end
             disp(norm(res)/npoints)
             if JFNKits==specs.JFNKitsMax
               disp('ran out of iterations, using fsolve')
                 opts=  optimset('display','off');
               %[uvals, fval,exitflag,output] = fsolve(@finite_diff_u_n1, uvals, opts);

%                if exitflag<0
%                    disp('Error: Solution not found')
%                    return;
%                end
               %disp('Oh no! Could not find a suiable w using linesearch')
               %break;
                sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
             end
             if norm(res)/npoints>0.1
                 for i = 1:1000
                     [Jac, Rhs]   = linearised_mom_op(sol, params, specs);
                     linearRes    = Jac*uvals - Rhs;
                     du           = Jac\(-linearRes);
                     [w, res, fl] = linesearch_w(sol, params, specs, du, linearRes);
                     if fl == 1
                         %uvals = uvals + w*du;
                         uvals = uvals + du;
                     else
                         disp('could not reduce residual, using fsolve 2!')
                         opts=  optimset('display','off');
                         %[uvals, fval,exitflag,output] = fsolve(@finite_diff_u_n1, uvals, opts);
                         uvals = uvals + du/32;
                         %                      if exitflag<0
                         %                          disp('Error: Solution not found')
                         %                          return;
                         %                      end
                         sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
                         %disp('Oh no! Could not find a suiable w using linesearch')
                         break;
                     end
                     sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];
                 end
             end


        case 'BISICLES'
                disp('Error: Velocity solver BISICLES does not exist yet')
            return;
    end

    sol = [Hvals; uvals; Nvals; phivals; hvals; helvals; h_cavvals; Lval];

%-------------------------------------------------------------------------%
%                   Residual functions
%-------------------------------------------------------------------------%

%Residual of momentum equation with n=1 (i.e. a Newtonian fluid)
    function res = finite_diff_u_n1(w)
        uvals = w(1:npoints);

        % Momentum equation
        mom_source = mutilde*max(Nvals(2:end-1).*(uvals(2:end-1)./(uvals(2:end-1) + lambdatilde(2:end-1).*Nvals(2:end-1))), 0) + beta*uvals(2:end-1)+ Hvals(2:end-1).*S_grad(2:end-1)/Lval;
        % boundary condition u=0 at z = 0
        k0    = uvals(1);
        % Centred difference discretisation of momentum equation
        k2nm1 = 4*epsilon*(H_midpoints(2:end).*(uvals(3:end) - uvals(2:end-1)) - H_midpoints(1:end-1).*(uvals(2:end-1) - uvals(1:end-2)))/dz^2 -Lval^2* mom_source;
        % Stress boundary condition at z = 1
        kn    = (Hvals(end)^2 -rhotilde_o*b_end^2*(b_end<0) )/8/Hvals(end)*dz + epsilon*(uvals(end-1) - uvals(end));

        % Residual vector
        res = [k0; k2nm1; kn];


    end

% Residual of momentum equation
    function res = finite_diff_u(w)
        uvals = w(1:npoints);

        
        % Momentum equation
        mom_source = mutilde*max(Nvals(2:end-1).*(uvals(2:end-1)./(uvals(2:end-1) + lambdatilde(2:end-1).*Nvals(2:end-1).^(n))).^(1/n), 0) + beta*uvals(2:end-1)+ Hvals(2:end-1).*S_grad(2:end-1)/Lval;
        % boundary condition u=0 at z = 0
        k0    = uvals(1);
        % Centred difference discretisation of momentum equation (regularised for small velocity gradients)
        k2nm1 = 4*epsilon*(H_midpoints(2:end).*(sqrt((uvals(3:end) - uvals(2:end-1)).^2 + regparam^2).^(1/params.n -1)).*(uvals(3:end) - uvals(2:end-1)) - H_midpoints(1:end-1).*(sqrt((uvals(3:end) - uvals(2:end-1)).^2  + regparam^2).^(1/params.n -1)).*(uvals(2:end-1) - uvals(1:end-2)))/dz^(1 + 1/n) -Lval^(1 + 1/n)* mom_source;
        % Stress boundary condition at z = 1
        kn    = (Hvals(end)^2 -rhotilde_o*b_end^2*(b_end<0) )/8/Hvals(end)*dz^(1/n)*Lval^(1/n) + epsilon*sqrt((uvals(end-1) - uvals(end))^2 + regparam^2)^(1/params.n -1)*(uvals(end-1) - uvals(end));

        
        %plot(4*epsilon*(H_midpoints(2:end).*(uvals(3:end) - uvals(2:end-1)) - H_midpoints(1:end-1).*(uvals(2:end-1) - uvals(1:end-2)))/dz^2)
        % Residual vector
        res = real([k0; k2nm1; kn]);


    end



end

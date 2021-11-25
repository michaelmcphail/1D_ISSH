classdef solverSpecs
  properties
     z_max          = 1;   % Domain size (z)
     npoints        = 300; % Number of z grid points
     t_max          = 1;   % Maximum time
     tpoints        = 100; % Number of time points (to export), assumes export time is uniform
     x_extended_max = 5;   % Length of spatial domain for which prescribed fields are known
     nx_points      = 1000; % Number of points on which to define prescribed fields
     t_current      = 0;    %current time

     regparam = 0.001;     % Viscosity regularisation parameter
     slope    = -0.5;      % Basal slope (eventually, b will be a function called in driver)
     artificialDiffusitivity = 0.1; %An artificial diffusivity term introduced into hcav evolution equation

     sol_0    = [];  %Initial guess

     %Calving mechanism {'constant', 'topographyDependent', 'flotationBased'}
     calvingType = 'constant'
     
     % Subglacial system selector (constitutive law for hel) {'poroelastic', 'exponential', 'logorithmic'}
     helConstitutive = 'poroelastic'
     
     %Surface melt function type
     msType          = 'periodicTanh'
     msPeriod        = 1;
     msSharpness     = 20;
     msAmplitude     = 0.5;
     
     periodicity_tol = 0.000001;
     
     

     %Chosen H solver, options are {'forwardsEuler', 'backwardsEuler', 'finiteVolumeUpwind', 'finiteVolumeFaceCentred'}
     HSolver  = 'forwardsEuler'

     %Chosen velocity solver, options are {'fsolveCentredDiff', 'JFNK', 'BISICLES'}
     usolver  = 'fsolveCentredDiff'

     %Velocity solver options
     JFNKtol               = 0.00001; % L2 tolerance defining convergence of JFNK solver
     JFNKitsMax            = 1000;    % Maximum number of iterations of JFNK solver before giving up
     linesearchAttemptsMax = 5;      % Maximum number of times to halve w (linesearch parameter) before giving up;

     %Continuation parameter, measure of how stationary L is
     dLTol   = 0.0005;
     dSolTol = 0.0005;

     %Chosen subglacial solver, options are {'centredDiff', 'JFNK'}
     subglacialSolver  = 'centredDiff'


     % Case: advance both glacier and subglacial hydrology?
     advanceGlacier    = 1; % Advance Glacier (H and u)? (1 = yes, 0 = no)
     advanceSubglacial = 1; % Advance Subglacial hydrology? (1 = yes, 0 = no)

     %Use real dimensionless parameters?
     realisticParams   = 0; %Use parameters in dimParams to calculate nonDimParams? (1 = yes, 0 = no)


     %Saving data
     saveData    = 1;  % Save data? (1 = yes, 0 = no)
     fieldsfn    = "./data/fields_tmax_n1_1.csv"; % File to save H and u
     tvalsfn     = './data/tvals_tmax_n1_1.csv'; %File to save tvalues


     % Plotting options
     plotDimensional = 0;  % Plot dimensional(0) or dimensionless (1) values
     increments      = 10; % Number of increments between plotted values
     figno           = 1;  % Figure number
     showplot        = 1;  % Show plot? (1 = yes, 0 = no)
     save            = 0;  % Save plot? (1 = yes, 0 = no)
     fname           = "all_fields_fixedN_n1.pdf";
  end
end

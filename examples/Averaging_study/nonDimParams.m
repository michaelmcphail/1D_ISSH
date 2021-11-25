classdef nonDimParams
   properties
       mbtilde  = 6*10^(-3); % Basal melt rate

       % Exponents present
       n        = 3;         % Glen's law exponent
       alpha    = 3;         % Hydraulic conductivity exponent
       gamma    = 0.5;       % h_el exponent

       %Momentum equation parameters
       epsilon     = 0.1;         % Dimensionless extensional viscosity coeff
       mutilde     = 1;           % Dimensionless basal friction coeff.
       beta        = 0.001;       % Dimensionless side friction coeff.
       lambdatilde = 4*10^(-3);   % Dimensionless ~roughness.
       rhotilde_w  = 1.1;         % Ratio of water to ice density.
       rhotilde_o  = 1.1;         % Ratio of ocean to fresh ice density.

       %subglacial parameters
       Ktilde      = 1;     % Dimensionless hydraulic conductivity
       sigmatilde  = 0.001; % Dimensionless rate of  change of englacial storage change
       delta       = 0.001; % Dimensionless rate of sheet height change
       nu          = 1;     % Ratio of Nscale to phiscale
       omega_r     = 1;     % ~Dimensionless rate of cavitation
       Atilde      = 1;     % Dimensionless rate of viscous
       hrtilde     = 1;     % Dimensionless basal roughness
       hctilde     = 1;     % Dimensionless coefficient of h_el
       mathcalV    = 0.5;   % Dimensionless strength of channel melt


       %Parameters related to calving
       C                     = 0.1; % Constant calving
       gamma_calving         = 2;   % Power law depending on depth
       thickness_sensitivity = 1;   % Sensitivity of ice to being above flotation value

       %Prescribed functions. Note: If not prescribed as functions of x
       %(hence z), say data instead, then we need to map the data to z
       %position.
       b_xvals   = linspace(0, 5, 1000); %x for topography
       a_xvals   = linspace(0, 5, 1000); %x for accumulation
       ms_xvals  = linspace(0, 5, 1000); %x for surface melt

       b  = [0.1 - 0.5* linspace(0, 0.5-0.005, 100).^2, 0.1 + 0.5/4 - 0.5*linspace(0.5, 5, 900)]; % data values for topography
       a  = 1*ones(1000, 1)'; % data values for accumulation
       ms = 1.05*[(1 - 4*(0.5 - linspace(0, 0.5, 100)).^2), ones(1, 900)] - 6*10^(-3);  % data values for surface melting
       ms_subOnly = 1.05*[(1 - 4*(0.5 - linspace(0, 0.5, 100)).^2), ones(1, 900)] - 6*10^(-3);  % data values for for varying subglacial independently of SMB

       %Surface melt parameters
       r      = 19.9;        % significance of summer
       rs     = 21.7;        % Height sensitivity
       tspr   = 3.6*10^(-4); % dimensionless time before spring
       taut   = 6.5*10^(-4); % dimensionless time before Autumn
       deltat = 5.7*10^(-5); % dimensionless time between seasons
       period = 9.8*10^(-4); % dimensionless period of surface melt osscilation
       
       %hel Constitutive law parameters
       h0     = 0;     % poroelastic sheet value for N=0
       h1     = 1;     % logarithm coefficient for hel
       N0     = 1;     % dimensionless reference pressure
       deltaN = 0.5;   % sensitivity of hel to changes in effective pressure
       
       %Metrics for assessing accuracy
       cumulative_mass_loss = [];
       sea_level_height     = [];

   end

end

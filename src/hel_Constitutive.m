function [helOut, fprimeN] = hel_Constitutive(sol, params, specs)
% HEL_CONSTITUTIVE return the poroelastic component of sheet thickness
%    sol: full solution vector at current time
%    params: contains the dimensionless parameters
%    specs: contains the sovler specifications. 

%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%

%identify the type of constitutive law to use
helConstitutiveLawType = specs.helConstitutive;

npoints = specs.npoints;             % Number of z points

% value of solution at $t$
Hvals     = sol(1:npoints);
uvals     = sol(npoints+1:2*npoints);
Nvals     = sol(2*npoints+1:3*npoints);
phivals   = sol(3*npoints+1:4*npoints);
hvals     = sol(4*npoints+1:5*npoints);
helvals   = sol(5*npoints+1:6*npoints);
h_cavvals = sol(6*npoints+1:7*npoints);
Lval      = sol(end);


%temporary values until params is augmented to include (1 for comparison, 0
%for paper).
h0     = params.h0;
% h0     = 0;
h1     = params.h1;
N0     = params.N0;
deltaN = params.deltaN;


%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%

switch helConstitutiveLawType
    case 'exponential'
        helOut = h0*exp(-Nvals/N0);    
    case 'logorithmic'    
        helOut = h0 + h1*log(Nvals./(N0 + deltaN*Nvals)); %TODO: h1 sign
    case 'none'
        helOut = 0.55*helvals.^0;
end

switch helConstitutiveLawType    
    case 'exponential'
        fprimeN =  -h0*exp(-Nvals/N0)/N0;       
    case 'logorithmic'
        fprimeN = -h1*(N0 + deltaN*Nvals).*(-deltaN*Nvals./(N0 + deltaN*Nvals).^2 + 1./(N0 + deltaN.*Nvals))./Nvals;
    case 'none'
        fprimeN = 0*helvals;
end




end

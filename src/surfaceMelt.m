function [msOut] = surfaceMelt(sol, params, specs, t, msOriginal)
% SURFACEMELT return a surface melt profile
%    sol: full solution vector at current time.
%    params: contains the dimensionless parameters.
%    specs: contains the solver specifications.
%    t: current time.
%    msOriginal: the original melt rate (sometimes useful for reference).

%-------------------------------------------------------------------------%
%                   Initialise
%-------------------------------------------------------------------------%


%Specifications
npoints = specs.npoints;             % Number of z points
z_values    = linspace(0, specs.z_max, npoints); %grid on which solution is defined

% value of solution at $t$
Hvals     = sol(1:npoints);
uvals     = sol(npoints+1:2*npoints);
Nvals     = sol(2*npoints+1:3*npoints);
phivals   = sol(3*npoints+1:4*npoints);
hvals     = sol(4*npoints+1:5*npoints);
helvals   = sol(5*npoints+1:6*npoints);
h_cavvals = sol(6*npoints+1:7*npoints);
Lval      = sol(end);

%Temporary parameter values

if isprop(params, 'r')
    r      = params.r;
    rs     = params.rs;
    deltat = params.deltat;
    tspr   = params.tspr;
    taut   = params.taut;
    period = params.period;
end

if isprop(params, 'xmean')
    %For now C is the amplitude
    xmean  = params.xmean;
end



%-------------------------------------------------------------------------%
%                   return ms based on specified type
%-------------------------------------------------------------------------%

switch specs.msType
    case 'constant'
        msOut = msOriginal;
        
    case 'periodicTanh'
        amp        = specs.msAmplitude;
        period     = specs.msPeriod;
        sharpness  = specs.msSharpness;
        scalingFac = 1 - amp*(1 +  tanh(sharpness*(mod(t,period)/period  - 0.25)))/2 + amp*(1 + tanh(sharpness*(mod(t,period)/period - 0.75)))/2;
        
        msOut      = scalingFac*msOriginal;
        
    case 'elevationDependent'
        %Get H at x values Fill rest of vector, find S at x values by adding b. 
        % Find ms at x points contained within ice sheet, set rest to final
        % value of ms (or zero?), return this (Note assumption here that b
        % and ms are defined at the same x values
        xvalsForHinterp = params.ms_xvals(params.ms_xvals<Lval);
        Hvals_atx       = interp1(z_values, Hvals, xvalsForHinterp/Lval);
        Svals_atx       = Hvals_atx + params.b(1:sum(params.ms_xvals<Lval));
        
        msOutShort      = max(0, r*(0.5*(tanh((mod(t, period) - tspr)/deltat) - tanh((mod(t, period) - taut)/deltat))) - rs*Svals_atx);
        msOut           = [msOutShort, msOutShort(end)*ones(1, length(msOriginal) - length(msOutShort))];
        
    case 'msSweep'
        % sinusoidal forcing around base meltrate
        xvalsForHinterp = params.ms_xvals(params.ms_xvals<Lval);
        Hvals_atx       = interp1(z_values, Hvals, xvalsForHinterp/Lval);
        Svals_atx       = Hvals_atx + params.b(1:sum(params.ms_xvals<Lval));
        
        msOutShort      = max(0, min(r - rs*Svals_atx, r))*(1 - specs.msAmplitude*sin(2*pi*t/period));
        msOut           = [msOutShort, msOutShort(end)*ones(1, length(msOriginal) - length(msOutShort))];
        
        
    case 'msSweepNonConserved'
        %Sinusoidal forcing with moving ablation point
        xvalsForHinterp = params.ms_xvals(params.ms_xvals<Lval);
        Hvals_atx       = interp1(z_values, Hvals, xvalsForHinterp/Lval);
        Svals_atx       = Hvals_atx + params.b(1:sum(params.ms_xvals<Lval));
        
        msOutShort      = max(0, min(r*(1 - specs.msAmplitude*sin(2*pi*t/period)) - rs*Svals_atx, r));
        msOut           = [msOutShort, msOutShort(end)*ones(1, length(msOriginal) - length(msOutShort))];
        
    case 'suddenMelt'
        %Sudden on--off switch for melt
        xvalsForHinterp = params.ms_xvals(params.ms_xvals<Lval);
        Hvals_atx       = interp1(z_values, Hvals, xvalsForHinterp/Lval);
        Svals_atx       = Hvals_atx + params.b(1:sum(params.ms_xvals<Lval));
        
        modt            = mod(t, period) /period;
        sharpness       = 30;
        modulationFac   = specs.msAmplitude*(-tanh(sharpness*modt)+ tanh(sharpness*(modt - 0.5)) - tanh(sharpness*(modt - 1)));

        

        msOutShort      = max(0, min(r - rs*Svals_atx, r))*(1 - modulationFac);
        msOut           = [msOutShort, msOutShort(end)*ones(1, length(msOriginal) - length(msOutShort))];        
        
    case 'pulse'
        % A single pulse of melt
        xvalsForHinterp = params.ms_xvals(params.ms_xvals<Lval);
        Hvals_atx       = interp1(z_values, Hvals, xvalsForHinterp/Lval);
        Svals_atx       = Hvals_atx + params.b(1:sum(params.ms_xvals<Lval));
        
        xmean           = params.C;
        twidth          = 0.1/60;
        tmean           = 0.1;
        xwidth          = Lval* 0.0005;
        pulse           = specs.msAmplitude.*exp(-(t - tmean).^2/twidth).*exp(-(xvalsForHinterp - xmean).^2/xwidth);
        
        msOutShort      = max(0, min(r - rs*Svals_atx, r)) + pulse;
        msOut           = [msOutShort, msOutShort(end)*ones(1, length(msOriginal) - length(msOutShort))];
        
        
end

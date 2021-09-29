function [rp,op,fp,strp,k,klabels,desc] = ...
  get_parameters(k,klabels,varargin)

  %Parameters used from Reynolds
  rp.nu1 = 0.08; % source of inflammatory cells; unit M*h^-1
  descs.nu1 = '% source of inflammatory cells; unit M*h^-1';
  rp.nu2 = 0.12; % decay of inflammatory cells; unit h^-1
  descs.nu2 = 'decay of inflammatory cells; unit h^-1';
  rp.muM = 0.05; % decay rate of inflammatory cells; unit h^-1
  descs.muM = 'decay rate of inflammatory cells; unit h^-1';
  rp.a1 = 0.04; %maximum production rate of anti-inf
  descs.a1 = 'nu3; maximum production rate of anti-inflammatory response; unit A*h^-1';
  alt_name.a1 = 'nu3';
  rp.muA = .05; %decay rate of the anti-inflammatory cells; unit h^-1 (Reynolds thought it was too high (Their value was .1)
  descs.muA = 'decay rate of the anti-inflammatory cells; unit h^-1 (Reynolds thought it was too high (Their value was .1))';
  rp.kM=.01; %acivation of inf by inf (value from reynolds)
  descs.kM = 'activation of inf by inf (value from reynolds); unit 1*(Mh)^-1';
  %  Adjusted from Reynolds so that kB*c1 = 0.1 = k_{np} from Reynolds'
  %  paper
  rp.kB=.1; %activation of inf by bacteria (value from reynolds)
  descs.kB = 'activation of inf by bacteria; unit 1*(10^6 cells/cc*h)^-1 (value from reynolds)';
  rp.kepsilon=.02; %activation of inf by damage (value from reynolds)
  descs.kepsilon = 'activation of inf by damage; unit 1*(eps*h)^-1 (value from reynolds)';
  rp.kAMEps=48; %relative effectivness of activated phagocytes and damaged tissue inducing production of the anti-inflammatory mediator
  descs.kAMEps = 'k4; relative effectivness of activated phagocytes and damaged tissue inducing production of the anti-inflammatory mediator; unit M/eps';
  alt_name.kAMEps = 'k4';
  rp.kBl=.6; %rate at which the non-specific local response eliminates pathogen (/L-units/h)
  descs.kBl = 'k2; rate at which the non-specific local response eliminates pathogen; unit 1/M/h';
  alt_name.kBl = 'k2';
  rp.sl=.005; %source of non-specific local response (L-units/h)
  descs.sl = 'source of non-specific local response; unit M/h';
  rp.mul=.002; %Decat rate for the non-specific local response (/h)
  descs.mul = 'Decay rate for the non-specific local response; unit 1/h';
  rp.klB=.01; %rate at which the non-specific local response is exhausted by pathogen (/B-units/h)
  descs.klB = 'k3; rate at which the non-specific local response is exhausted by pathogen; unit 1/B/h';
  alt_name.klB = 'k3';
  rp.sA=.0125; %source of anti-inflammatory (A-units /h)
  descs.sA = 'source of anti-inflammatory; unit A/h';

  %Other parameters
  %  Alli's
  op.tau=24; % time scale for epithelium repair; unit h
  descs.tau = 'time scale for epithelium repair; unit h';
  op.c1=1; %virulence***No longer in use/set to 1***
  descs.c1 = 'virulence of pathogen; strength with which bacteria instigate an immune response; unitless';
  op.epsilon0=0; %baseline damage
  descs.epsilon0 = 'baseline damage; eps';
  op.T =1; % Threshold value for damage; unit cells2g^-1*h^-1
  descs.T = 'Threshold value for damage; unit eps/h';
  op.k1=.5; %growth rate for bacteria
  descs.k1 = 'growth rate for bacteria; unit 1/h';
  %  Purely for finding steady state estimates...if integration takes
  %  longer than this amount of time; integration is halted and a warning
  %  issued that we may not have yet reached steady state.
  op.maxintegrationtime = 10;
  descs.maxintegrationtime = 'max cpu time that we allow the ode solver to run for; unit cpu s';
  op.minintegrationtime = 0.01;
  descs.minintegrationtime = 'min cpu time that we allow the ode solver to run for; unit cpu s';
  op.gamma = 0;0.5;
  descs.gamma = 'used for estimating kA with gamma = 0 meaning 1/(1+kA*Amax) = 25% while gamma = 1 means (1/(1+kA*Amax))/(1/(1+gamma*kA*Amax)) = 25%; unitless';
  if op.gamma > 0
    error('Nonzero gamma not implemented (see solver_for_kA)');
  end
  %  Corresponds to 75% inhibition
  op.inhi_perc = 0.75;
  descs.inhi_perc = '% inhibition that we desire (see gamma); unitless';

%   op.tau=0.1; % time scale for epithelium repair; unit h
%   op.c1=0.5; %virulence
%   op.epsilon0=0; %baseline damage
%   op.T =1; % Threshold value for damage; unit cells2g^-1*h^-1,
%   op.k1=.5; %growth rate for bacteria

%   %  Amy's
%   op.tau = 4;
%   op.c1 = 0.5;
%   op.epsilon0 = 0;
%   op.T = 1;
%   op.k1 = 0.5;

%  %  Amy's ode file
%  rp.kepsilon = 0.05;
%  

  %  Set default values for these unknown parameters
%   fp.Dexp = 100;
%   fp.Dmax = 100;
  fp.kA=7.1834182161733215466;1/0.28;2.5960;
  descs.kA = 'inhibition rate of the anti-inflammatory response; unit 1/A';
  fp.f=15.3782;
  descs.f = 'maximum rate of damage produced by the pro-inflammatory response; unit eps/M/h';
  fp.k5=1.8;%1.6236;
  descs.k5 = 'rate at which pro-inflammatory response consumes pathogen; unit 1/M/h';
  fp.B_infy=145.0376;
  descs.B_infy = 'pathogen carrying capacity';
  %  Not currently used
%   fp.k_De=0.097329351856;
  fp.Bsource_in = 1;
  descs.Bsource_in = 'initial bacterial load injected into peritoneum; 10^6 bacteria';
  %  Jared introduced 12/3/2018 to check out rescaling
  %  Modified 3/19/2019
  fp.tscale = 1;
  descs.tscale = 'rescales time…not sure if this is currently fully operational; unit; 1/h';
  %  Jared introduced 5/1/2021.  Experiments were supposedly allowed to
  %  "breathe" for 24 h before measurements started (at time 0h)
  fp.td = 0;
  descs.td = 'time delay-how long we wait before "starting" the experiment';
  %  6/18/2019-We separate these up into two different rescaling factors,
  %  Bmf and Bsf that convert bacterial values from the simulation into
  %  "measured values" (Bmf-associated with the measured values in the
  %  excel sheet) and bacterial values from the simulation into "source
  %  values" (Bsf-associated with Bsource/initial values of bacteria in the
  %  clot).  Mathematically speaking; Bsf = Bmf = 1 would be ideal.
  %  Physiologically speaking; Bsf and Bmf should be determined by
  %  consulting with our collaborators.  Having had difficulties with that,
  %  we will go ahead and let our optimizer find optimal values for these
  %  parameters that allow us to best fit the data.
  fp.Bmf = 1;
  descs.Bmf = '1/Bs; for rescaling bacteria with Bs*B = 10^6 bacteria/cc; unit cc';
  alt_name.Bmf = '1/Bs';
  %  This is equal to Bmf/Bsf.
  fp.Bmfsf = 1;
  descs.Bmfsf = 'not used but once was something like the measured units (from blood measurements)/Bsource units (from estimated peritoneal injection #s; units measured bacteria/source bacteria';
  %  6/18/2019-Note that changing Mscale; Ascale; and escale the way that
  %  we have proposed (varold -> varscale*varnew) should not change the
  %  bacterial levels 
%   fp.Mscale = 1;
%   fp.Ascale = 1;
%   fp.escale = 1;
  %  Setting clot_capacity to a huge value runs the code as normal.
  %  Setting it to a smaller value in some programs will cause B to be
  %  initially set to Bsource_in-clot_capacity.
  fp.clot_capacity = 2000;
  descs.clot_capacity = 'not currently functional (clot capacity is 2000e8; huge); Bsource units';
  fp.nu4 = 400;
  descs.nu4 = 'Half-saturation of anti-inflammatory response; unit M/eps';
%   fp.kB = 0.2;
  fp.kD1 = 1.5227e4;
  descs.kD1 = 'kD; decay of bacterial population in the clot (released into body); unit 1/h';
  alt_name.kD1 = 'kD';
  fp.kD2 = 0;
  descs.kD2 = 'used for more complicated dosing functions/models including bacteria in the clot as an additional modeled variable (5 odes instead of 4)';
  fp.kD3 = 0;
  descs.kD3 = 'used for more complicated dosing functions/models including bacteria in the clot as an additional modeled variable (5 odes instead of 4)';
  fp.kD4 = 0;
  descs.kD4 = 'used for more complicated dosing functions where D(t) is a shifted Gaussian';
  fp.use_log_Bc = false;
  descs.use_log_Bc = 'allows logistic growth of bacteria in a clot; currently unused';
  fp.limit_flux = false;
  descs.limit_flux = 'whether or not to limit bacteria flux at the clot/blood interface in a certain way (only can exit; can never enter clot?)';
  fp.effective_carrying_capacity = false;
  descs.effective_carrying_capacity = 'whether to use the effective carrying capacity B_infyxB_mf in place of the "carrying capacity" B_infy';
  
  %  Info for classification.  "a" = absolute; "r" = relative
  fp.Binfasep = 1e-6; %  Above this at t = inf; we have "septic"
  descs.Binfasep = 'criterion for classifying as healthy/aseptic/septic';
  fp.dinfasep = 1;  %  Above this at t = inf; we have "septic"|"aseptic"
  descs.dinfasep = 'criterion for classifying as healthy/aseptic/septic';
  fp.Bfinrsep = 2;  %  Above this (0.1 = 10%) sim still "septic"
  descs.Bfinrsep = 'criterion for classifying as healthy/aseptic/septic';
  %  This is the target maximum fractional value that aseptic simulations
  %  are at when they "finish" as determined by fp.fint.  Assigning a value
  %  of 1/eps makes it so that this parameter effectively never gets used.
  fp.Bfinasepmax = 1/eps;  
  descs.Bfinasepmax = 'criterion for classifying as healthy/aseptic/septic';
  %  This is the target minimum fractional value that septic simulations
  %  are at when they "finish" as determined by fp.fint.  Assigning a value
  %  of eps makes is so that this parameter effectively never gets used.
  fp.Bfinsepmin = eps;
  descs.Bfinsepmin = 'criterion for classifying as healthy/aseptic/septic';
  fp.fint = 200;  %  The "finite time" at which we check finite criteria
  descs.fint = 'final integration time…time at which we try to classify cases (vs infinity)';
  fp.sM=0; %source of pro-inflammatory (M-units /h)
  descs.sM = 'source of pro-inflammatory; unit M/h';
  
  %  Used to choose which dosing function to use
  %  Originally meant for "string" parameters...later on I shoved logical
  %  parameters in there too.
  strp.dosing_function = 'simple_exponential';
  strp.solve_for_kA = false;
  strp.check_kA = false;
  %  Added to consider subsets
  
  %Rewrites the values of the parameters in k
  for kc = 1:numel(klabels)
    if isfield(fp,klabels{kc})
      fp.(klabels{kc}) = k(kc);
    elseif isfield(op,klabels{kc})
      op.(klabels{kc}) = k(kc);
    elseif isfield(rp,klabels{kc})
      rp.(klabels{kc}) = k(kc);
    end
  end
  
  %  Varargin stuff put here
  strucnames = {'rp','op','fp','strp'};
  %  Loop through varargin cell
  for vac = 1:2:numel(varargin)
    %  Loop through structure names
    for sc = 1:numel(strucnames)
      strucname = strucnames{sc};
      eval(['tmp = ',strucname,';']);
      if isfield(tmp,varargin{vac})
        eval([strucname,'.',varargin{vac},' = varargin{vac+1};']);
      end
    end
  end
  
  %  Calculate alternate dosing function quantities (these quantities
  %  depend on kD1-4 so if these change...calculations need to be redone)
  if (fp.clot_capacity < 2000) && isequal(strp.dosing_function,...
      'shifting_gaussian')
    warning(['Alternate dosing function not ',...
      'set up to work with clot_capacity']);
  end
  if (fp.Bsource_in > 0) && isequal(strp.dosing_function,...
      'shifting_gaussian')
    fp.tD = fp.kD1*fp.Bsource_in^(-fp.kD2);
    descs.tD = 'center of the shifting Gaussian dose function';
    fp.Dmax = max((-fp.kD3)*fp.tD+fp.kD4,1e-3);
    descs.Dmax = 'maximum associated with the shifting Gaussian dose function';
    fp.Dexp = fzero(@(Dexp) fp.Bsource_in-fp.Dmax*...
      integral(@(t) exp(-Dexp.*(t-fp.tD).^2),0,inf),[10*eps,1e10]);
    descs.Dexp = 'coefficient in the exponent for the shifting Gaussian dose function';
    % fp.Dexp = 1;
  else
    fp.tD = 0; 
    descs.tD = 'center of the shifting Gaussian dose function';
    fp.Dmax = 0;
    descs.Dmax = 'maximum associated with the shifting Gaussian dose function';
    fp.Dexp = 0;
    descs.Dexp = 'coefficient in the exponent for the shifting Gaussian dose function';
  end
  
  %  Glean all klabels and the corresponding default values
  fnr = fieldnames(rp); fno = fieldnames(op); fnf = fieldnames(fp);
  klabels = {fnr{:},fno{:},fnf{:}};
  if numel(klabels) ~= numel(unique(klabels))
    error(['Duplicate definitions (probably) in Reynolds (rp),',...
      ' Other (op); and Free (fp) parameters']);
  end
  charinds = [];
  knew = [];
  klabelsnew = {};
  for kc = 1:numel(klabels)
    if isfield(rp,klabels{kc})
      tmp = rp.(klabels{kc});
    elseif isfield(op,klabels{kc})
      tmp = op.(klabels{kc});
    elseif isfield(fp,klabels{kc})
      tmp = fp.(klabels{kc});
    end
    if ~ischar(tmp)
      knew(end+1) = tmp;
      klabelsnew{end+1} = klabels{kc};
    end
  end
  k = knew;
  klabels = klabelsnew;
  
  %  Solving for kA added 8/19/2019
  [fp2,k,klabels] = solver_for_kA(rp,op,fp,strp,k,klabels,1940);
  fp.Bmax = fp2.Bmax;
  if strp.solve_for_kA
    fp.kA = fp2.kA;
  end
  
  desc.descs = descs;
  desc.alt_name = alt_name;

%   fprintf('kA = %30.20g; k5 = %30.20g.\n',fp.kA,fp.k5);
end
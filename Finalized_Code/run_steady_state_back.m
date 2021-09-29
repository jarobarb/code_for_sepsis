% This function runs the simulation and returns the steady state values of y
% For example, to run with A0=0.5 and nu=0.1: run_steady_state('A0',0.5,'nu1',0.1)
function [ystar,t,y] = run_steady_state(varargin)

  %  Default parameter values
  odetol = 1e-2;
  mult_param_vals = 1;
  mult_param_name = 'dummie';

  %Parameters
  [k,klabels] = define_default_ks;
  for vac = 1:2:numel(varargin)
    switch varargin{vac}
      case 'k', k = varargin{vac+1};
      case 'klabels', klabels = varargin{vac+1};
      case 'mult_param_vals', mult_param_vals = varargin{vac+1};
      case 'mult_param_name', mult_param_name = varargin{vac+1};
    end
  end
%   [rp,op,fp,strp] = get_parameters(k,klabels,varargin{:});
  if isempty(mult_param_vals)
    [rp,op,fp,strp] = get_parameters(k,klabels,varargin{:});
    rp(1) = rp; op(1) = op; fp(1) = fp; strp(1) = strp;
  else
    for mpc = 1:numel(mult_param_vals)
      [rp(mpc),op(mpc),fp(mpc),strp(mpc)] = ...
        get_parameters(k,klabels,varargin{:},...
        mult_param_name,mult_param_vals(mpc));
    end
  end

  % some error handling
  if ~mod(numel(varargin),2)==0
    error('inputs must be parameter-value pairs such as ''k1'',0.5')
  end

%   % Rewrites over selected parameter values
%   for p=1:2:numel(varargin)
%     evalc([varargin{p};]); % returns an error if you enetered a parameter that does not exist
%     eval([varargin{p},'=varargin{p+1};']);
%   end

  %  5/17/2019: A0 = 0.25 eliminated in favor of A0 = Asteadystate = sA/muA
%   if isempty(mult_param_vals)
%     B0=0.0; M0=0.0; A0=rp.sA/rp.muA; Eps0=0.0;
%     yhealthy = [B0 M0 A0 Eps0];
%     B0 = fp.Bmfsf*max(fp.Bsource_in-fp.clot_capacity,0);
%     %   fp.Bsource_in = fp.Bsource_in-B0;
%     y0=[B0,M0,A0,Eps0];
%     for vac = 1:2:numel(varargin)
%       if isequal(varargin{vac},'y0'), y0 = varargin{vac+1}; end
%       if isequal(varargin{vac},'B0'), y0(1) = varargin{vac+1}; end
%       if isequal(varargin{vac},'M0'), y0(2) = varargin{vac+1}; end
%     end
%     y=y0;
%   else
    y = []; yh = [];
    warning('Multiple parameter values not enabled for B0, M0, A0, Eps0');
    for mpc = 1:numel(mult_param_vals)
      B0=0.0; M0=0.0; A0=rp(mpc).sA/rp(mpc).muA; Eps0=0.0;
      yhealthy = [B0 M0 A0 Eps0];
      B0 = fp(mpc).Bmfsf*max(fp(mpc).Bsource_in-fp(mpc).clot_capacity,0);
      %   fp.Bsource_in = fp.Bsource_in-B0;
      y0=[B0,M0,A0,Eps0];
      for vac = 1:2:numel(varargin)
        if isequal(varargin{vac},'y0'), y0 = varargin{vac+1}; end
        if isequal(varargin{vac},'B0'), y0(1) = varargin{vac+1}; end
        if isequal(varargin{vac},'M0'), y0(2) = varargin{vac+1}; end
      end
      for vac = 1:2:numel(varargin)
        if isequal(varargin{vac},'y0'), y0 = varargin{vac+1}; end
      end
      y=[y,y0]; yh = [yh,yhealthy];
    end
    y0 = y;
%   end
  
  use_old = false;

  Prmc.tmax=1/eps;%1000000; % max time to let the simulation run
  Prmc.NNmax=100;Prmc.tmax/10; % sets th enumber of intervals over which you check to see if at steady state. I set NNMax so that each interval is always 50 time steps
  Prmc.fulltime=false; % set to true to make the simulation run the fulltime
  tolvec=1e-8.*ones(size(y)); % tolerance for the standard deviation of each population at steady state
%   tolvec=5e-4.*ones(1,4); % tolerance for the standard deviation of each population at steady state

  if use_old
    for NN = 1:Prmc.NNmax

      tvec = Prmc.tmax*[NN-1,NN]./Prmc.NNmax; % tvec=[start time end time] for current interval

      % Stops if any population is less than zero
      if any (y(end,:)<-1e-10)
        warning('some population is negative');
        disp('current y-values [B,M,Eps,A]:');
        fprintf('%8f ',y(end,:));
        fprintf('\nsimulation time: %f\n',tvec(1));
        for p=1:2:numel(varargin)
          fprintf(strcat(varargin{p},'=%f\n'),varargin{p+1});
        end
      end

      %     [~,y] = ode15s(dydt,tvec,y(end,:),...
      %       odeset('InitialStep',0.1,'NonNegative',ones(size(y(end,:))),...
      %       'MaxStep',1000,'RelTol',0.1*min(tolvec),'AbsTol',0.1*min(tolvec)));

      %  Adaptive time step
      [t,y] = ode15s(@dydt, tvec, y(end,:), ...
        odeset('InitialStep',0.1,'NonNegative',ones(size(y(end,:))),...
        'MaxStep',1000,'RelTol',odetol*min(tolvec),...
        'AbsTol',odetol*tolvec),rp,op,fp,strp);

%       [t,y,te,ye,ie] = ode15s(@dydt, tvec, y(end,:), ...
%         odeset('InitialStep',0.1,'NonNegative',ones(size(y(end,:))),...
%         'MaxStep',1000,'RelTol',0.01*min(tolvec),...
%         'AbsTol',0.01*min(tolvec),'Jacobian',@myjac,...
%         'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,tolvec)),...
%         rp,op,fp,strp);

      avgyall = mean(y,1); % returns row vector with the average value of each variable over the last tvec
      stdyall = mean(abs(y-avgyall),1); % returns row vector with std of each variable over the last tvec

      stdyallquant = sum(stdyall < tolvec); % returns the number of variables whose std satisfies tolvec

      % If all the variables satisfy tolvec, then it's at steady state!
      if (stdyallquant == numel(avgyall)) && (~Prmc.fulltime)
        ystar=y(end,:);
        break
      elseif (NN==Prmc.NNmax)
        warning('did not get to steady state for tmax=%f',Prmc.tmax);
        fprintf('\ninitial y0 = \t');
        fprintf('%f ',y0);
        fprintf('\ncurrent y = \t');
        fprintf('%8f ',y(end,:));
        fprintf('\nstdyall:\t\t');
        fprintf('%f ',stdyall);
        disp(' ');
        for p=1:2:numel(varargin)
          fprintf(strcat(varargin{p},'=%f\n'),varargin{p+1});
        end
        ystar=5.*ones(4,1); % makes ystar high so mapping will decrease par2
        %error('did not get to steady state for tmax=%f',Prmc.tmax');
      end
    end
  else

    %  Alternate route via JOB, 7/1/2019
    %  Check to make sure this isn't B_sourcein = 0 (should always be a
    %  health state)
    if (all([fp.Bsource_in] == 0)) && (isequal(y0,yh))
      ystar = y0; y = ystar; t = 0;
    else
      tstart = tic;
      opts = odeset('InitialStep',odetol,...'NonNegative',0*ones(size(y(end,:))),...
        'MaxStep',10000,'RelTol',odetol*min(tolvec),...
        'AbsTol',odetol*min(tolvec),'Jacobian',@mymultjac,...
        'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
        odetol*tolvec,tstart));
      %  For set time step (comment/uncomment as need be)
%       dtset = 10;
%       opts = odeset('InitialStep',dtset,...'NonNegative',0*ones(size(y(end,:))),...
%         'MaxStep',dtset,'RelTol',1/eps,...
%         'AbsTol',1/eps,'Jacobian',@myjac,...
%         'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
%         0.01*tolvec,tstart));
      [t,y,te,ye,ie] = ode15s(@dydt,[0,Prmc.tmax],y0,opts,rp,op,fp,strp);
%       fprintf('dts: ');
%       fprintf('%g ',unique(diff(t)));
%       fprintf('\n');
      ystar = y(end,:);
      if ie == 2 %t(end) == Prmc.tmax
%         warning('Did not get to steady state for tmax=%f',Prmc.tmax);
        warning(['Did not get to steady state after integrating for ',...
          '%g seconds!'],op.maxintegrationtime);
        for p=1:2:numel(varargin)
          if isfloat(varargin{p+1})
            fprintf(strcat(varargin{p},'=%f\n'),varargin{p+1});
          end
        end
        ystar = 5*ones(4,1);
      end
    end
  end
end
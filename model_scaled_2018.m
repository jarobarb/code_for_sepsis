function [sepsis_exp,sol,rp,op,fp,strp,k,klabels] = ...
  model_scaled_2018(varargin)

  %  k-vector of values of the free parameters
  %  klabels-cell of names of the free parameters
  %  by passing in both these arguments we can relatively easily pick and
  %  choose which free parameter values we want to explore
  %  time:  times that the bacteria in blood are calculated to compare
  %  with experimental measurement times
  
  %  model_scaled_2018 is made for only 1 bacteria at a time, hence
  %  Bsource_ins (plural) doesn't make sense
  for vac = 1:2:numel(varargin)
    if isequal('Bsource_ins',varargin{vac})
      varargin{vac} = 'Bsource_in';
      if numel(varargin{vac+1}) > 1
        error('Multiple Bsources for model_scaled_2018');
      end
    end
  end
  
  %  Default parameter values
  overalltol = 1e-2;
  
  plot_for_paper = true;

  [k,klabels] = define_default_ks;
  Bsource_in = 1;  
  tspan = [0 200];
  time = [];
  outlier_matrix = [];

  %y(1)= Bacteria;
  %y(2)= Pro-inflammatory population;
  %y(3)= Anti-inflammatory population;
  %y(4)= Damage;
  
  %  Set a default value for plotting (i.e. if false, don't plot)
  plot_solution = true;

  %  Whether or not classify the corresponding run (as septic, aseptic, or
  %  healthy...see my_finder and run_steady_state...those two codes should
  %  agree with here).
  classify = false;
  
  %  For checking out the individual terms of the bacteria eqn
  plot_b_terms = false;

  %  Rewrite over the above default value (i.e. if true for plot_solution,
  %  go ahead and plot the solution now)
  for vac = 1:2:numel(varargin)
    if exist(varargin{vac}) == 1
      eval([varargin{vac},' = varargin{vac+1};']);
    end
  end
  
  %"Graphical" parameters
  %  Line color
  gp.lc = 'b';
  %  Line style
  gp.ls = '-';
  %  Matlab lies about how it renders figures.  If you say I want a 5x6
  %  figure, it somehow converts that into pixels, inaccurately, and uses
  %  that.  This factor theoretically fixes that issue and is specific to
  %  each computer.
  gp.cf = 1.167;%6/5.25;
  %  Line width
  gp.lw = 1;
  %  Figure number
  gp.fn = 1;
  for vac = 1:2:numel(varargin)
    if isfield(gp,varargin{vac})
      gp.(varargin{vac}) = varargin{vac+1};
    end
  end
  %  Font size
  gp.fs = 12*gp.cf;

  %Parameters
  [rp,op,fp,strp,k,klabels] = get_parameters(k,klabels,varargin{:});
  
  %  Data
  if ~exist('exd','var')
    [exd] = get_data(outlier_matrix,'td',fp.td);
  end
  
  %  Initial conditions
%   Bsource_ind = ismember(Bsource_in,exd.Bsource);
%   if isempty(Bsource_ind)
%     error('You picked a Bsource for which there is no experimental data');
%   else
    fp.Bsource_in = Bsource_in;
%   end
  % Initial conditions
  %  Adjusted:  A0 = sA/(muA*Ascale) (not 0.5)...  
  y0 = [max(fp.Bsource_in-fp.clot_capacity,0); 0; rp.sA/rp.muA; 0; min(fp.clot_capacity,fp.Bsource_in)]; % y0(4)=0.1 from page 3 of the first paper%%%%%%

  for vac = 1:2:numel(varargin)
    if isequal(varargin{vac},'y0'), y0 = varargin{vac+1}; end
  end
  
  %System of differential equations
%   [t,y] = ode15s(@dydt,[tspan(1),tspan(end)],...
%     y0,odeset('NonNegative',true),rp,op,fp,strp);
  tstart = tic;
  tolvec=1e-8.*ones(1,5);
  %  *** Matches run_steady_state ***
  %  The adaptive time stepper makes sure |e(i)| <=
  %  max(RelTol*abs(y(i)),AbsTol(i)).  If we pick RelTol super small
  %  (small enough so that RelTol*abs(y(i)) is always < AbsTol(i)) and also
  %  pick AbsTol(i) big enough in the middle (we try 1e10), then only the
  %  first entry (when considering multiple parameter values) and last
  %  entry should control the step size.  This will, ideally, make the time
  %  stepping the exact same for the middle parameter values so long as the
  %  first and last parameter values stay the same.
  abstolvec=1e-8.*ones(size(tolvec)); 
  % tolerance for the standard deviation of each population at steady state
  sstolvec = 1e-8.*ones(size(tolvec));
  reltol=1e-8;eps;
  if classify
%     opts = odeset('InitialStep',0.1,'NonNegative',ones(size(y0)),...
%       'MaxStep',10000,'RelTol',0.01*min(tolvec),...
%       'AbsTol',0.01*min(tolvec),'Jacobian',@myjac,...
%       'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
%       0.01*tolvec,tstart));
    opts = odeset('InitialStep',overalltol,...'NonNegative',0*ones(size(y(end,:))),...
        'MaxStep',10000,'RelTol',overalltol*reltol,...
        'AbsTol',overalltol*min(abstolvec),'Jacobian',@mymultjac,...
        'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
        overalltol*sstolvec,tstart));
%     opts = odeset('InitialStep',0.1,'NonNegative',ones(size(y0)),...
%       'MaxStep',10000,'RelTol',0.01*min(tolvec),...
%       'AbsTol',0.01*min(tolvec),...
%       'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
%       0.01*tolvec,tstart));
%     opts = odeset('InitialStep',overalltol,...'NonNegative',0*ones(size(y(end,:))),...
%         'MaxStep',10000,'RelTol',overalltol*reltol,...
%         'AbsTol',overalltol*min(abstolvec));%,'Jacobian',@mymultjac);
    odes_sol = ode15s(@dydt,[0,1/eps],y0,opts,rp,op,fp,strp);
    sol.classification = my_classifier(odes_sol);
  else
%     odes_sol = ode15s(@dydt, [tspan(1),tspan(end)],...
%       y0,odeset('InitialStep',0.1,'NonNegative',ones(size(y0)),...
%       'MaxStep',1000,'RelTol',0.01*min(tolvec),...
%       'AbsTol',0.01*min(tolvec),'Jacobian',@myjac),rp, op, fp, strp);
    odes_sol = ode15s(@dydt, [tspan(1),tspan(end)],...
      y0,odeset('InitialStep',0.1,'NonNegative',ones(size(y0)),...
      'MaxStep',1000,'RelTol',0.01*min(tolvec),...
      'AbsTol',0.01*min(tolvec)),rp, op, fp, strp);
  end
  
  t = odes_sol.x; y = odes_sol.y';
  if ~isempty(time)
    kok = time < t(end);
    knotok = time > t(end);
    if any(kok)
      tmp = deval(odes_sol,time(kok));%interp1(t,y(:,1),time);
    else
      tmp = zeros(4,0);
    end
    sepsis_exp = y(1,end)*ones(size(time));
    sepsis_exp(kok) = tmp(1,:);
  else
    sepsis_exp = [];
  end
  sol.solution_time = t;
  sol.bacteria = y(:,1);
  sol.pro_inf = y(:,2);
  sol.anti_inf = y(:,3);
  sol.damage = y(:,4);
  sol.clot = y(:,5);
  sol.eq3a = y(:,2)+rp.kAMEps*y(:,4);
  sol.aprod = rp.a1*(sol.eq3a)./(1+fp.kA*y(:,3))./...
    (fp.nu4+sol.eq3a);
  est_flux = fp.kD1*sol.clot-fp.kD2*sol.bacteria;
  flux_ind = (~fp.limit_flux) | ((sol.clot < fp.clot_capacity) | ...
    (est_flux > 0));
  sol.bflux = est_flux.*flux_ind;
  sol.struc = odes_sol;
  sol.struc.rp = rp; sol.struc.op = op; sol.struc.fp = fp;
  sol.struc.strp = strp;

  if plot_solution
    figure(gp.fn);  set(gcf,'Units','inches','Position',[1,1,6,4]);
    
    if plot_for_paper, subplot(2,2,1); else, subplot(2,3,1); end
    plot(t,y(:,1),[gp.lc,gp.ls],'Linewidth',gp.lw*gp.cf);
    if max(y(:,1)) > 4000000 %40
      yup = 0.4;
      ylow = 40;
      ylim([0,90]);
      breakinfo = breakyaxis([0.4,40]);
      breakinfo.lowAxes.YTick = [0,0.3];
      breakinfo.highAxes.YTick = [50,90];
%       lowpos = breakinfo.lowAxes.Position;
%       highpos = breakinfo.highAxes.Position;
%       breakinfo.breakAxes.YLim = [0,1];
%       breakinfo.lowAxes.Position = lowpos+[0 0 0 0.1];
%       breakinfo.highAxes.Position = highpos+[0 0.1 0 -0.1]
%       keyboard
    end
    
    fonnam = 'Times New Roman';
    fonnam = 'Arial';
    
    hold on;
    ylabel('Bacterial Levels (10^6 bacteria)','FontSize',gp.fs,'FontName',fonnam)
    xlabel('Time (hrs)','FontSize',gp.fs,'FontName',fonnam)
    set(gca,'FontSize',gp.fs,'FontName',fonnam,'LineWidth',gp.lw);

    if plot_for_paper, subplot(2,2,2); else, subplot(2,3,2); end
    plot(t,y(:,2),[gp.lc,gp.ls],'Linewidth',gp.lw*gp.cf);
    hold on;
    ylabel('Pro-inflammatory Levels','FontSize',gp.fs,'FontName',fonnam)
    xlabel('Time (hrs)','FontSize',gp.fs,'FontName',fonnam);
    set(gca,'FontSize',gp.fs,'FontName',fonnam,...
        'LineWidth',gp.lw);

    if plot_for_paper, subplot(2,2,3); else, subplot(2,3,3); end
    plot(t,y(:,3),[gp.lc,gp.ls],'Linewidth',gp.lw*gp.cf);
    hold on;
    ylabel('Anti-inflammatory Levels','FontSize',gp.fs,'FontName',fonnam)
    xlabel('Time (hrs)','FontSize',gp.fs,'FontName',fonnam);
    set(gca,'FontSize',gp.fs,'FontName',fonnam,...
        'LineWidth',gp.lw);

    if plot_for_paper, subplot(2,2,4); else, subplot(2,3,4); end
    plot(t,y(:,4),[gp.lc,gp.ls],'Linewidth',gp.lw*gp.cf);
    hold on;
    ylabel('Damage Levels','FontSize',gp.fs,'FontName',fonnam)
    xlabel('Time (hrs)','FontSize',gp.fs,'FontName',fonnam);
    set(gca,'FontSize',gp.fs,'FontName',fonnam,...
        'LineWidth',gp.lw);
      
    if ~plot_for_paper
      subplot(2,3,5)
      plot(t,y(:,5),[gp.lc,gp.ls],'Linewidth',gp.lw*gp.cf);
      hold on;
      ylabel('Clot bacteria','FontSize',gp.fs,'FontName',fonnam)
      xlabel('Time (hrs)','FontSize',gp.fs,'FontName',fonnam);
      set(gca,'FontSize',gp.fs,'FontName',fonnam,...
        'LineWidth',gp.lw);
    end

    %     figure
    %     hold on
    %     plot(t,y(:,1),'b','Linewidth',2);
    %     ylabel('Bacteria')
    %     xlabel('time')
    %     hold on
    %     plot(mtime,conc, 'r*');
    %     legend('Estimated Bacteria Level','Actual Bacteria Level')


    %{
  figure
  plot(Bsource,D1(Bsource), 'Linewidth',2);
  ylabel('D1');
  xlabel('Bsource');

  figure
  plot(Bsource,D2(Bsource), 'Linewidth',2);
  ylabel('D2');
  xlabel('Bsource');
    %}

    %{
  figure
  plot(t,D(fp.Bsource_in, t), 'Linewidth',2);
  ylabel('D');
  xlabel('time');
    %}

    %{

  plot(t,y(:,1) ,'Linewidth',2);
  hold on
  ylabel('Bacteria')
  xlabel('time')
  plot(mtime,conc, '*');
  hold on
  legend('128: Estimated','128: Actual', '248: Estimated','248: Actual', '505: Estimated','505: Actual', '1940: Estimated','1940: Actual')

  figure
  plot(time1,actual_conc1,'b')
  ylabel('bacteria in blood')
  xlabel('time')
  hold on
  plot(time2,actual_conc2,'g')
  plot(time3,actual_conc3,'r')
  plot(time4,actual_conc4,'k')
  legend('Bsource=1.28e8','Bsource=2.48e8','Bsource=5.05e8','Bsource=1.94e9')
    %}

    %{
   figure
   plot(t,fM, 'Linewidth',2);
   ylabel('')
   xlabel('time')
   hold on
   plot(t, ones(size(t))*fp.T,'Linewidth',2);
   legend('fM','T')
   hold off
    %}

  end
  
  if plot_b_terms
    B = sol.bacteria; M = sol.pro_inf; A = sol.anti_inf;
    carrying_capacity = fp.effective_carrying_capacity*fp.B_infy+...
      (~fp.effective_carrying_capacity)*fp.B_infy*fp.Bmf;
    sol.Blog = op.k1.*B.*(1-B/carrying_capacity);
    sol.loci = rp.kBl*rp.sl*B./(rp.mul+rp.klB*(B/fp.Bmf));
    sol.gloi = fp.k5*B.*M./(1+fp.kA*A);
    if plot_solution
      figure(10);
      plot(sol.solution_time,op.k1.*B.*(1-B/carrying_capacity),...
        sol.solution_time,rp.kBl*rp.sl*B./(rp.mul+rp.klB*(B/fp.Bmf)),...
        sol.solution_time,fp.k5*B.*M./(1+fp.kA*A));
      hold on
    end

  end

end
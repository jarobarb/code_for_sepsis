% function [bacteria_blood_exp, sol] = sepsis(k,klabels,time,Bavg,fp.Bsource_in,varargin)
function plot_all_points_2018(varargin)

  %  k-vector of values of the free parameters
  %  klabels-cell of names of the free parameters
  %  by passing in both these arguments we can relatively easily pick and
  %  choose which free parameter values we want to explore
  %  time:  times that the bacteria in blood are calculated to compare
  %  with experimental measurement times
  
  overalltol = 1e-2;

  [k,klabels] = define_default_ks;
  outlier_matrix = [];
  plot_what = 'ave';
  
  %  Bsource vector
  col = 'gbrc';
  
  %  Rewrite over the above default value (i.e. if true for plot_solution, go
  %  ahead and plot the solution now)
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  %y(1)= Bacteria;
  %y(2)= Pro-inflammatory population;
  %y(3)= Anti-inflammatory population;
  %y(4)= Damage;

  %  Parameters
  [rp,op,fp,strp] = get_parameters(k,klabels,varargin{:});
  
  %  Initial conditions
  %  Adjusted:  A0 = sA/(muA*Ascale) (not 0.5)...  
  y0 = [max(fp.Bsource_in-fp.clot_capacity,0); 0; rp.sA/rp.muA; 0; min(fp.clot_capacity,fp.Bsource_in)]; % y0(4)=0.1 from page 3 of the first paper%%%%%%

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
  
  %  Data
  [exd] = get_data(outlier_matrix,'td',fp.td);
  
  opts = odeset('InitialStep',overalltol,...'NonNegative',0*ones(size(y(end,:))),...
    'MaxStep',10000,'RelTol',overalltol*reltol,...
    'AbsTol',overalltol*min(abstolvec),'Jacobian',@mymultjac,...
    'Events',@(t,y,rp,op,fp,strp) myevent(t,y,rp,op,fp,strp,...
    overalltol*sstolvec,tstart));

  for Bloadc=1:numel(exd.Bsource)
    fp.Bsource_in = exd.Bsource(Bloadc);
    tfinal = max(exd.t.allo{Bloadc});
    t = linspace(0,tfinal);
    tolvec=1e-8.*ones(1,4);
    % Diff eq solving
    y0(1) = fp.Bmfsf*max(fp.Bsource_in-fp.clot_capacity,0);
    y0(end) = min(fp.clot_capacity,fp.Bsource_in);
    
%     y0(1) = max(fp.Bsource_in-fp.clot_capacity,0);
%     fp.Bsource_in = fp.Bsource_in-y0(1);
    [t,y] = ode15s(@dydt,t,y0,opts,rp,op,fp,strp);
    
    sol.solution_time = t;
    sol.bacteria = y(:,1);
    sol.pro_inf = y(:,2);
    sol.anti_inf = y(:,3);
    sol.damage = y(:,4);

    % Plots the model and the data points
    %  Matlab lies about how it renders figures.  If you say I want a 5x6
    %  figure, it somehow converts that into pixels, inaccurately, and uses
    %  that.  This factor theoretically fixes that issue and is specific to
    %  each computer.
    compfact = 1.167;
    lw = compfact*1;
    fons = compfact*12;
    fonnam = 'Times New Roman';
    plot(t,y(:,1) ,'color', col(Bloadc), 'Linewidth', lw)
    hold on
    set(gca,'FontSize',fons,'LineWidth',lw)
    ylabel('Bacterial Levels (10^6 bacteria)','fontsize',fons,...
        'FontName',fonnam);
    xlabel('Time (hours)','fontsize',fons,'FontName',fonnam);
    switch plot_what
      case 'ave'
        p1 = plot(exd.t.ave{Bloadc},exd.B.ave{Bloadc},[col(Bloadc),'x']);
        p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
      case 'all'
        for expc = 1:numel(exd.t.all{Bloadc})
          for tc = 1:numel(exd.t.all{Bloadc}{expc})
            if ~exd.outlier{Bloadc}{expc}(tc)
              p1 = plot(exd.t.all{Bloadc}{expc}(tc),...
                exd.B.all{Bloadc}{expc}(tc),[col(Bloadc),'.'],...
                'MarkerSize',15);
              p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            else
              p1 = plot(exd.t.all{Bloadc}{expc}(tc),...
                exd.B.all{Bloadc}{expc}(tc),[col(Bloadc),'o']);
              p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
          end
        end
    end
    for legc = 1:numel(exd.Bsource)
      legcell{legc} = ['B_{source} = ',num2str(exd.Bsource(legc)),...
        '\times{10^6} bacteria'];
    end
    legend(legcell{:},'Location','best');
    legend('boxoff');
    
  end

end
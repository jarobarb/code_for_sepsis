function [cat,tc,maxdam,sol,dampts] = mortality_estimate(varargin)

  %  cat-"category" that the simulation falls into.  Referring to mortality
  %  data jpg, red is the "dies in < 24 hr" category", blue-"dies in < 48
  %  hr", yellow-"dies in < 72 hr", black-"dies in < 96 hr"
 
  %  varargin-a big long list of parameter value pairs likely to include k,
  %  klabels (see "define_default_ks"/"get_parameters", for instance) as
  %  well as other parameter pairs (see mortality_script_for_Julia for
  %  typical call)
  
  crit_val = 25;
  %  For mortality data, the last time check is 96 hrs
  maxt = 96;
  fignum = 1;
  %  Line color
  lc = 'b';
  
  for vac = 1:2:numel(varargin)
    if exist(varargin{vac})
      eval([varargin{vac},' = varargin{vac+1};']);
    end
  end
  
  [sepsis_exp,sol,rp,op,fp,strp] = model_scaled_2018(varargin{:},...
    'tspan',[0,maxt],'fn',fignum);
  
  maxdam = max(sol.damage);
  
  if maxdam < crit_val
    tc = inf;
  else
    %  Note, implicit (interp1) assumption of damage solution's
    %  monotonicity...e.g. curve only crosses crit_val once
    %  extrap should not be necessary, we just always include it for safe
    %  measure
    %  The first few damage entries can, in fact, be repeated 0s, which
    %  interp1 doesn't like.  Get rid of them.  Unless crit_val is super
    %  small (between 0 and the first non-zero entry, which is usually
    %  O(odesolvertol)), there should be no harm done.
    [c,ia,ic] = unique(sol.damage);
    tc = interp1(sol.damage(ia),sol.solution_time(ia)',crit_val,...
      'linear','extrap');
  end
  
  if tc < 24
    cat = 1;
  elseif tc < 48
    cat = 2;
  elseif tc < 72
    cat = 3;
  elseif tc < 96
    cat = 4;
  else
    %  Can't really tell that this is a special case from experimental data
    %  so we make it into cat 4 again.
    cat = 4;
  end
  
  for timec = [24,48,72,96]
    dampts.(['d',num2str(timec)]) = interp1(sol.solution_time,...
      sol.damage,timec,'linear','extrap');
  end

end
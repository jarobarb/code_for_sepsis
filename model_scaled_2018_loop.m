function [val,sepsis_exp,sol,rp,op,fp,strp,kcell,klabelcell] = model_scaled_2018_loop(...
  Bsource,t,exd,k,klabels,varargin)

  parcomp = false;
  %  For classification, we use the same classification found in
  %  "my_finder" and the same integration as found in "run_steady_state"
  classify_vec = false(size(exd.Bsource));

  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  if exist('fixedk','var')
    k = [k(:);fixedk(:)]';
    klabels = {klabels{:},fixedklabels{:}};
  end
  
  if isinparallel | (~parcomp)
    for Bloadc = 1:numel(exd.Bsource)
      [sepsis_exp{Bloadc},sol{Bloadc},rp{Bloadc},op{Bloadc},fp{Bloadc},...
        strp{Bloadc},kcell{Bloadc},klabelcell{Bloadc}] = model_scaled_2018('k',k,...
        'klabels',klabels,'time',exd.t.ave{Bloadc},...
        'Bsource_in',exd.Bsource(Bloadc),...
        'classify',classify_vec(Bloadc),varargin{:});
    end
  else
    parfor Bloadc = 1:numel(exd.Bsource)
      [sepsis_exp{Bloadc},sol{Bloadc},rp{Bloadc},op{Bloadc},fp{Bloadc},...
        strp{Bloadc},kcell{Bloadc},klabelcell{Bloadc}] = model_scaled_2018('k',k,...
        'klabels',klabels,'time',exd.t.ave{Bloadc},...
        'Bsource_in',exd.Bsource(Bloadc),...
        'classify',classify_vec(Bloadc),varargin{:});
    end  
  end
  
  if isempty(Bsource)
    val = [];
  else
    %  The two vectors should be the same size
    val = zeros(size(Bsource));
    for Btc = 1:numel(Bsource)
      Bvec = exd.Bsource;
      for Bloadc = 1:numel(exd.Bsource)
        if t(Btc) <= max(sol{Bloadc}.solution_time)
          tmp = deval(sol{Bloadc}.struc,t(Btc));
          Bvec(Bloadc) = tmp(1);
        else
          Bvec(Bloadc) = sol{Bloadc}.bacteria(end);
        end
      end
      val(Btc) = interp1(exd.Bsource,Bvec,Bsource(Btc),'linear','extrap');
    end
  end

end
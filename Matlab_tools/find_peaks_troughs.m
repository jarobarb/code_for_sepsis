function [pk_vals,pk_locs,tr_vals,tr_locs] = find_peaks_troughs(x,varargin)

  %  smaller peaks within minpkdist of larger peaks will be destroyed
  minpkdist = 0;
  minpkheight = [];
  %  smaller peaks within 2*minpkdist of larger peaks will be destroyed if
  %  they are (after shifting of x) less than adjpkheight*large peak height
  adjpkheight = 0;
  
  for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  xmean = mean(x);
  x = x-xmean;
  
  if ~isempty(minpkheight), pk_lim_inds = find(x>minpkheight); end
  
  trend = sign(diff(x));
  flat_inds = find(trend == 0);
  
  N = length(trend);
  
  if ~isempty(flat_inds)
    if flat_inds(end) < N-1
      trend(flat_inds(end)) = trend(flat_inds(end)+1);
    end
    
    for fic = length(flat_inds)-1:-1:1
      trend(flat_inds(fic)) = trend(flat_inds(fic)+1);
    end
  end
  
  ind_maxs = find(diff(trend)==-2)+1;
  
  if ~isempty(minpkheight)
    pk_locs = intersect(pk_lim_inds,ind_maxs);
  else
    pk_locs = ind_maxs;
  end
  pk_vals = x(pk_locs);
  
  if minpkdist > 0
    
    [pk_vals,tmp_inds] = sort(pk_vals,'descend');
    pk_locs = pk_locs(tmp_inds);
    
    pk_ind_del = false(size(pk_locs));
    
    for pc = 1:length(pk_locs)
      if ~pk_ind_del(pc)
        pk_ind_del = pk_ind_del | ...
          ( (pk_locs >= pk_locs(pc)-minpkdist) &...
          (pk_locs <= pk_locs(pc)+minpkdist) );
        pk_ind_del = pk_ind_del | ...
          ( (pk_locs >= pk_locs(pc)-2*minpkdist) &...
          (pk_locs <= pk_locs(pc)+2*minpkdist) &...
          (pk_vals <= adjpkheight*pk_vals(pc)) );
        pk_ind_del(pc) = false;
      end
    end
    pk_locs(pk_ind_del) = [];
    pk_vals(pk_ind_del) = [];
    
    [pk_locs,tmp_inds] = sort(pk_locs,'ascend');
    pk_vals = pk_vals(tmp_inds);
    
  end
  pk_vals = pk_vals+xmean;
  
  %  Now find associated minimums
  %  First find halfway points between peaks
  if length(pk_vals) > 1
    tr_vals(length(pk_vals)-1) = 0;
    for pc = 1:length(pk_vals)-1
      [tr_vals(pc),tmp_ind] = min(x(pk_locs(pc):pk_locs(pc+1)));
      tr_locs(pc) = pk_locs(pc)+tmp_ind-1;
    end
  else
    tr_vals = []; tr_locs = [];
  end
   
%   pk_mid_locs = round((pk_locs(2:end)+pk_locs(1:end-1))./2);
%   
%   [tr_vals(length(pk_locs)),tmp_ind] = min(x(pk_mid_locs(end):end));
%   tr_locs(length(pk_locs)) = tmp_ind+pk_mid_locs(end)-1;
%   [tr_vals(1),tmp_ind] = min(x(1:pk_mid_locs(1)));
%   tr_locs(1) = tmp_ind;
%   for pc = 1:length(pk_mid_locs)-1
%     [tr_vals(pc+1),tmp_ind] = min(x(pk_mid_locs(pc):pk_mid_locs(pc+1)));
%     tr_locs(pc+1) = pk_mid_locs(pc)+tmp_ind-1;
%   end

end
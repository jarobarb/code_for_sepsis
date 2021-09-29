function [xval] = par_bisec(fnct,xbnds,varargin)

  %  Inputs:
  %    fnct-function we want to find the zero of
  %    xbnds-leftmost and rightmost x that we believe there is a zero
  %      between
  %    varargin-list of parameter values-pairs
  %  
  %  Outputs:
  %    xval-estimated location of the root

  %  Initialize frequently used variables.
  xreltol = 1e-12;
  xabstol = 1e-12;
  freltol = 1e-12;
  fabstol = 1e-12;
  parcomp = true;
  showmons = true;

  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  if parcomp
    p = gcp();
    n = p.NumWorkers;
  else
    n = 1;
  end
  
  xs = linspace(xbnds(1),xbnds(2),n);
  fs = zeros(size(xs));

  if parcomp
    for isc = 1:numel(xs)
      my_jobs(isc) = parfeval(p,fnct,1,xs(isc));
    end
    for isc = 1:numel(xs)
      [job_index,value] = fetchNext(my_jobs);
      fs(job_index) = value;
      if showmons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
    end
  else
    for isc = 1:numel(xs)
      fs(isc) = feval(fnct,xs(isc));
      if showmons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
    end
  end
  %  Use this to estimate the typical size of f
  fscale = max(abs(fs));
  
  if max(fs)*min(fs) > 0
    error('Zero not bracketed!  Widen/change xbnds and try again.');
  end
  
  while (min(diff(xs)) > xabstol) & (min(abs(fs)) > fabstol) & ...
      (min(diff(xs)) > xreltol*max(mean(xs),xabstol)) & ...
      (min(abs(fs)) > freltol*fscale)
    
    indr = max(find(fs > 0,1),find(fs < 0,1));
    indl = indr-1;
    
    xs = linspace(xs(indl),xs(indr),n+2);
    tmp = zeros(size(xs));
    tmp(1) = fs(indl); tmp(2) = fs(indr);
    fs = xs; fs(1) = tmp(1); fs(end) = tmp(2);
    
    if parcomp
      for isc = 1:(numel(xs)-2)
        my_jobs(isc) = parfeval(p,fnct,1,xs(isc+1));
      end
      for isc = 1:(numel(xs)-2)
        [job_index,value] = fetchNext(my_jobs);
        fs(job_index+1) = value;
        if showmons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
      end
    else
      for isc = 1:(numel(xs)-2)
        fs(isc+1) = feval(fnct,xs(isc+1));
        if showmons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
      end
    end
    
    if showmons
      fprintf('Dx = %g, miny = %g.\n',...
        min(diff(xs)),min(abs(fs)));
    end

  end
  
  xval = interp1(fs,xs,0,'linear','extrap');

end
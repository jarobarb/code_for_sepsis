function [xlow,ylow] = par_opt1d(xl,xm,xr,fnct,varargin)

  %  Initialize frequently used variables.
  xtol = 1e-3;
  ftol = 1e-6;
  abstol = 1e-6;
  parcomp = true;
  show_mons = true;

  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  if parcomp
    p = gcp();
    n = p.NumWorkers;
    n2 = ceil(n/2);
  else
    n = 1;
    n2 = 1;
  end

  xreldiff = 2*xtol; freldiff = 2*ftol;
  xs = linspace(xl,xr,max(3,n-1));
  xs = union(xs,xm);
  ys = xs;

  if parcomp
    for isc = 1:numel(xs)
      my_jobs(isc) = parfeval(p,fnct,1,xs(isc));
    end
    for isc = 1:numel(xs)
      [job_index,value] = fetchNext(my_jobs);
      ys(job_index) = value;
      if show_mons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
    end
  else
    for isc = 1:numel(xs)
      ys(isc) = feval(fnct,xs(isc));
      if show_mons, fprintf('Job #%d/%d done\n',isc,numel(xs)); end
    end
  end
  ylowini = min(ys);

  while (xreldiff > xtol) & (freldiff > ftol)

    [~,ind] = min(ys);
    if (ind == 1) || (ind == numel(ys))
      num_of_nearest_pts = 2;
    else
      num_of_nearest_pts = 3;
    end
    inds = 1:numel(ys);
    [a,b] = mink((ind-inds).^2,num_of_nearest_pts); b = sort(b);
    xsold = xs(b); ysold = ys(b);

    if numel(b) == 2
      xs = linspace(min(xsold),max(xsold),n+2);
      ys = nan(size(xs));
      ys(1) = ysold(1); ys(end) = ysold(2);
    else
      xs = union(linspace(xsold(1),xsold(2),n2+2),...
        linspace(xsold(2),xsold(3),n2+2));
      ys = nan(size(xs));
      ys(1) = ysold(1); ys(end) = ysold(3); ys(n2+2) = ysold(2);
    end

    xreldiff = max(diff(xsold))./max(mean(abs(xsold)),abstol);
    freldiff = max(diff(ysold))./max(mean(abs(ysold)),abstol);

    if parcomp
      clear('my_jobs');
      inds = find(isnan(ys));
      for ic = 1:numel(inds)
        xc = inds(ic);
        my_jobs(ic) = parfeval(p,fnct,1,xs(xc));
      end
      for ic = 1:numel(inds)
        yc = inds(job_index);
        [job_index,value] = fetchNext(my_jobs);
        yc = inds(job_index);
        ys(yc) = value;
        if show_mons, fprintf('Job #%d/%d done\n',ic,numel(xs)); end
      end
    else
      for xc = 1:numel(xs)
        if isnan(ys(xc))
          ys(xc) = feval(fnct,xs(xc));
        end
        if show_mons, fprintf('Job #%d/%d done \n',xc,numel(xs)); end
      end

    end
    
    [ylow,ind] = min(ys);
    xlow = xs(ind);
    if (ylow > ylowini) | any(isnan(ys))
      keyboard
    end
    
    if show_mons, fprintf('Dx = %g, Dy = %g, xlow = %g, ylow = %30.20g\n',...
        xreldiff,freldiff,xlow,ylow); end

  end

end
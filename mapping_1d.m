function [low,high] = mapping_1d(varargin)

  parcomp = true;
  m1dguessmin = 0;
  m1dguessmax = 505;
  m1dtol = 1e-9;eps;%1e-12;
  myparamname = 'Bsource_in';
  
  for vac = 1:2:numel(varargin)
    if exist(varargin{vac}) == 1
      eval([varargin{vac},' = varargin{vac+1};']);
    end
  end
  
%   %  Converts measured bacterial values to simulated bacterial values
%   m1dguessmin = m1dguessmin/Bmfsf;
%   m1dguessmax = m1dguessmax/Bmfsf;
  
  %  If we are inside another parallelized process, don't try to
  %  parallelize this process, Matlab will fail if we try to
  if isinparallel, parcomp = false; end
  
  if parcomp
    p = gcp();
    nproc1 = max(p.NumWorkers);
  else
    nproc1 = 1;
  end
  nproc2 = max(1,floor(nproc1/2));
  
  low = m1dguessmin;
  high = m1dguessmax;
  guesses = linspace(m1dguessmin,m1dguessmax,max(nproc1,3));
  uvals = nan(size(guesses));
  mexit = false;
  
%   while (min(diff(guesses)) > m1dtol/Bmfsf) || (min(diff(guesses)) == 0)
  while ~mexit%(min(diff(guesses)) > m1dtol) || (min(diff(guesses)) == 0)
    
    if parcomp
      for gc = 1:numel(guesses)
        if isnan(uvals(gc))
          my_jobs(gc) = parfeval(p,@my_finder,1,varargin{:},myparamname,...
            guesses(gc));
        end
      end
      for gc = 1:numel(guesses)
        if any(isnan(uvals))
          [job_index,value] = fetchNext(my_jobs);
          uvals(job_index) = value;
        end
%         fprintf('Job #%d/%d done\n',job_index,numel(guesses));
      end
    else
      for gc = 1:numel(guesses)
        if isnan(uvals(gc))
          uvals(gc) = feval(@my_finder,varargin{:},myparamname,guesses(gc));
          if isinf(uvals(gc)), keyboard; end
        end
%         fprintf('Job #%d/%d done\n',gc,numel(guesses));
      end
    end
   
    inds = find(diff(uvals));
    if any(diff(uvals) < 0) %(numel(inds) > 2) | (uvals(end) < max(uvals))
      warning(['Non-monotonic behavior wrt Bsource_in direction!  ',...
        'Using narrowest aseptic window guesstimate.']);
      fprintf('Tolerances = %g %g\n',guesses(end/2)-guesses(1),...
        guesses(end)-guesses(end/2+1));
      low = max(guesses(1:end/2)); high = min(guesses((end/2+1):end));
%       low = m1dguessmin; high = m1dguessmax;
      return
    end
    if (uvals(end) == 1) & (numel(guesses) == 3)
      warning(['Highest Bsource not septic death!  ',...
        'Setting to septic death.']);
      uvals(end) = 2;
      inds = find(diff(uvals));
    end
    locdiffs = diff(guesses); locdiffs = locdiffs(inds);
    if max(locdiffs) < m1dtol, mexit = true; end
    switch numel(inds)
      case 1
        guesses = linspace(guesses(inds),guesses(inds+1),nproc1+2);
        uvalsold = uvals;
        uvals = nan(size(guesses));
        uvals(1) = uvalsold(inds); uvals(end) = uvalsold(inds+1);
      case 2
        guesses = [linspace(guesses(inds(1)),guesses(inds(1)+1),nproc2...
          +2),linspace(guesses(inds(2)),guesses(inds(2)+1),nproc2+2)];
        uvalsold = uvals;
        uvals = nan(size(guesses));
        uvals(1) = uvalsold(inds(1));
        uvals(nproc2+2) = uvalsold(inds(1)+1);
        uvals(nproc2+3) = uvalsold(inds(2));
        uvals(end) = uvalsold(inds(2)+1);
      %  All values of uvals should be the same
      case 0
        switch unique(uvals)
          case 0
            low = high; return
          case 2
            high = low; return
          otherwise
            keyboard
        end
      otherwise
        keyboard
%         error('strange behavior!');
    end
    
  end
  
  finalguesses = guesses(~isnan(uvals));
%   if numel(inds) == 2
%     low = Bmfsf*mean(finalguesses(1:2));
%     high = Bmfsf*mean(finalguesses(3:4));
%   else
%     warning('Two bifurcation values not found');
%     if max(uvals) == 1
%       low = Bmfsf*mean(finalguesses);
%     else
%       high = Bmfsf*mean(finalguesses);
%     end
% %     keyboard;
%   end
  if numel(inds) == 2
    low = mean(finalguesses(1:2));
    high = mean(finalguesses(3:4));
  else
    warning('Two bifurcation values not found');
    if max(uvals) == 1
      low = mean(finalguesses);
    elseif min(uvals) == 1
      high = mean(finalguesses);
    else
      low = mean(finalguesses);
      high = low;
    end
%     keyboard;
  end
  
end
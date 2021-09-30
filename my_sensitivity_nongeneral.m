function [dudk,uvals_norm,klabels,kmat,uvals,s_struc] = ...
  my_sensitivity(k,klabels,addlargs,varargin)

  %  This program can be used to estimate the sensitivity of utility
  %    values with respect to various parameter values.  After the
  %    parameter values to be tested have been determined, the program uses
  %    simple first order finite differences along each parameter direction
  %    (effectively evaluating the utility function at n+1 points where n =
  %    number of parameters to be tested) to estimate the sensitivities.
  %
  %  INPUTS
  %  k-numerical values of parameters that we wish to be tested.
  %    Sensitivity calculations will be done near this point in parameter
  %    space.
  %  klabels-the names of the parameters...necessary for many of the
  %    inflammation programs
  %  addlargs-any additional arguments that need to be handed into the
  %    utility function
  %
  %  OUTPUTS
  %  varargin-accepts parameter pair values including the "default"
  %    my_sensitivity program specific ones as well as parameter value
  %    pairs for the actual ode model that we plan to eventually run
  %  dudk-dutilityfunction/dparameter-vector of estimates of sensitivity of
  %    the utility function wrt the parameters (numel(dudk) = number of
  %    parameters tested)
  %  dutil-dutilityfunction by itself.  Since all parameters are
  %    (currently) bumped the exact same relative amount, dutil gives a
  %    quick estimate as to which parameters are most/least sensitive.
  %    Again, vector of size numel(dutil) = number of parameters tested.
  %  klabels-names of the parameters that were eventually tested
  %  kmat-the n+1 points in parameters space listed in matrix form
  %  uvals-the n+1 utility function values for each point tested in
  %    parameter space

  %%  Sensitivity analysis
  %  This depends on the location of the point in parameter space as well as
  %  the magnitude of the bumps taken to estimate the underlying
  %  derivatives/sensitivities
  
  %%  Reshape k and klabels so all goes smoothly downstairs
  [mk,nk] = size(k);
  k = k(:);
  klabels = {klabels{:}}';

  %%  (Re)define default values for this program
  %  Default relative bump to use when estimating sensitivities of the
  %  utility function.
  sp.frac = sqrt(eps);
  %  Default utility function for testing
  sp.util = @(k,klabels,addlargs) sum(k.^2);
%   sp.util = @(k,klabels,addlargs) optk_scaled_2018(k,klabels,addlargs{:});
  sp.parcomp = false;
  %  Which performance metric to use (see optk_scaled_2019 stats structure
  %  for possible options)
  sp.metric = 'curragreement';
  %  If > 0, this uses a random normal distribution with std given below
  %  for checking out sensitivity.  Otherwise, we use structured sampling.
  sp.std = sqrt(eps);
  sp.npts = numel(k)*3;
  %  Plots marginal distributions
  sp.plotmds = true;
  %  Try some model fits
  sp.modelfit = true;
  
  %  Redefine default values
  for vac = 1:2:numel(varargin)
    sp.(varargin{vac}) = varargin{vac+1};
  end
  
  %%  Run the set of points in parameter space
  if sp.std <= 0
    kmat = repmat(k,1,numel(k)+1);
    kmat(:,2:end) = kmat(:,2:end)+diag(k.*sp.frac);
%     kmat(:,end+[1:numel(k)]) = kmat(:,2:end)-2*diag(k.*sp.frac);
  else
    kmat = repmat(k,1,sp.npts);
    kmat(:,2:end) = kmat(:,2:end)+diag(k.*sp.std)*...
      randn(numel(k),sp.npts-1);
  end

  myit = 1;
  if sp.parcomp
    %  Use parallel computation
    p = gcp();
    for kc = 1:(size(kmat,2))
      if sp.std <= 0
        mydir = (-1)^(kc > (numel(k)+1)).*sp.frac*([1:numel(k)]'==...
          (mod(kc-1,numel(k)+1)+(kc > numel(k)+1)));
      else
        mydir = sp.std;
      end
      my_jobs(kc) = parfeval(p,@myf,3,sp.util,kmat(:,1),klabels,...
        addlargs,mydir);
    end    
    uvals = zeros(1,size(kmat,2));
    for kc = 1:size(kmat,2)
      [job_index,value1,struc1,sol1] = fetchNext(my_jobs);
%       uvals(job_index) = value1;
      strucs(job_index) = struc1;
      uvals(job_index) = struc1.(sp.metric);
      sols{job_index} = sol1;
      kmat(:,job_index) = struc1.k;
%       fprintf('Job #%d/%d done\n',job_index,size(kmat,2));
      fprintf('Job #%d/%d done\n',kc,size(kmat,2));
    end
  else
    %  Without parallel computation
    for kc = 1:size(kmat,2)
      if sp.std <= 0
        mydir = sp.frac*([1:numel(k)]'==(kc-1));
        if kc > 1, fprintf('%s\n',klabels{kc-1}); end
      else
        mydir = sp.std;
      end
      [uvals(kc),strucs(kc),sols{kc}] = ...
        myf(sp.util,kmat(:,1),klabels,addlargs,mydir);
      uvals(kc) = strucs(kc).(sp.metric);
      kmat(:,kc) = strucs(kc).k;
      if isinf(uvals(kc)), keyboard; end
      fprintf('Job #%d/%d done\n',kc,numel(k)+1);
    end
  end
  
  %%  Use the results from above to calculate sensitivities
  %  
%   if sp.std <= 0
% %     dutil = uvals(2:end)'-uvals(1);
%     for sc = 2:(numel(strucs)-numel(k))
%       if ((strucs(sc).badSSs > 5) && (strucs(sc+numel(k)).badSSs > 5))
%         keyboard
%       elseif (strucs(sc).badSSs > 5)
%         smmat_shif(sc-1,:) = strucs(1).simmeas-strucs(sc+numel(k)).simmeas;
%         dutil(sc-1) = uvals(1)-uvals(sc+numel(k));
%       elseif (strucs(sc+numel(k)).badSSs > 5)
%         smmat_shif(sc-1,:) = strucs(sc).simmeas-strucs(1).simmeas;
%         dutil(sc-1) = uvals(sc)-uvals(1);
%       else
%         smmat_shif(sc-1,:) = strucs(sc).simmeas-strucs(sc+numel(k)).simmeas;
%         dutil(sc-1) = uvals(sc)-uvals(sc+numel(k));
%       end
%     end
%     %  Relative sensitivity of each measurement to each parameter:  du/dk
%     %  \approx (u(k+k*delta)-u(k))/(k*delta) so that k/u*du/dk = ...
%     %  (u(k+k*delta)-u(k))/(u(k)*delta)
%     dsmmat = smmat_shif./(strucs(1).simmeas(:)'.*sp.frac);
%     s_struc.smsqr = sqrt(sum(dsmmat.^2,2)./size(dsmmat,2));
%     s_struc.dsmmat = dsmmat;
%     s_struc.smmat_shif;
%     %  For debugging:
%     %  close(figure(373)); figure(373); pcolor(smat); colorbar
%     dudk = dutil./(uvals(1).*sp.frac);
%   else
    penalized = false(size(uvals));
    for uc = 1:numel(uvals)
      if (strucs(uc).badSSs > 5)
        penalized(uc) = true;
      end
    end
    fprintf('Penalized runs = %g\n',sum(penalized));
    uvals = uvals(~penalized);
    kmat = kmat(:,~penalized);
    strucs = strucs(~penalized);
    for sc = 1:numel(strucs)
      smmat(sc,:) = strucs(sc).simmeas;
    end
    %  We shift to hopefully give better numerical results in general
    if sp.std <= 0
      mean_smmat = smmat(1,:);
      mean_kmat = kmat(:,1);
    else
      mean_smmat = mean(smmat);
      mean_kmat = mean(kmat,2);
    end
    smmat_norm = (smmat-mean_smmat)./mean_smmat;
    kmat_norm = [ones(1,sum(~penalized));...
      (kmat-mean_kmat)./max(mean_kmat,eps(mean_kmat))];
    %  kmat_shif'*coefs = smmat_shif comes from the linear equation
    %  c0j*1+c1j*k1i+c2j*k2i+...=mij.  Taking the derivative of this
    %  expression with respect to each parameter yields c1j, c2j, etc. as
    %  the estimated derivative values.
    if sp.std > 0
      coefs = pinv(kmat_norm')*(smmat_norm);
      dsmmat = coefs(2:end,:);
    else
%       kmat_norm_tmp = kmat_norm(2:end,:);
      dk_norm = diag(kmat_norm(2:end,2:end));
%       coefs = 0*pinv(kmat_norm')*(smmat_norm); coefs = coefs(2:end,:);
%       for kc = 1:numel(klabels)
%         coefs = coefs+pinv(kmat_norm_tmp(:,1+[0,kc+[0,numel(k)]])')*...
%           smmat_norm(1+[0,kc+[0,numel(k)]],:);
%       end
%       dsmmat = coefs;
      dsmmat = smmat_norm(2:end,:)./dk_norm;
    end
    
    %  We do the same for the "measurement" u, the cost function
    uvals_norm = (uvals(:)-mean(uvals))./mean(uvals);
    coefs = (kmat_norm')\(uvals(:)-mean(uvals));
    dudk = coefs(2:end);
    
    s_struc.dsmmat = dsmmat;
    s_struc.smsqr = sqrt(sum(dsmmat.^2,2)./size(dsmmat,2));
    s_struc.smmat_shif = smmat_norm;
    s_struc.kmat_shif = kmat_norm;
%   end
  
  %%  Reshape everything to match the dimensions of the original k
  k = reshape(k,mk,nk);
  klabels = reshape(klabels,mk,nk);
  dudk = reshape(dudk,mk,nk);
%   uvals_norm = reshape(uvals_norm,mk,nk);
  
  %%  Plot results via marginal distributions
%   for kc = 1:numel(klabels)
%     eval([klabels{kc},' = kmat(kc,:)'';']);
%   end
%   uvalsc = uvals';
%   eval(['tbl = table(',sprintf('%s,',klabels{:}),'uvalsc);']);
%   
%   for kc = 1:numel(klabels)
% %     mdl1 = stepwiselm(tbl,'constant','ResponseVar','uvalsc');
%     eval(['tbl1 = table(',sprintf('%s,',klabels{kc}),'uvalsc);']);
%     mdl1 = fitlm(tbl1);
%     fprintf('%10s %g\n',klabels{kc},mdl1.Coefficients.pValue(end));
% %     if sp.plotmds
% %       tmpinteger = ceil(sqrt(numel(k))./sqrt(2));
% %       for j = 1:numel(k)
% %         subplot(tmpinteger,2*tmpinteger,j);
% %         scatter(kmat(j,:),uvals,'.');
% %       end
% %     end
%   end
  
end

function [myuval,mystruc,mysols] = myf(f,kmean,klabels,addlargs,std)

  myc = 0;
  if numel(std) > 1
    k = kmean.*(1+std.*ones(size(kmean)));
  else
    k = kmean.*(1+std.*randn(size(kmean)));
  end
  [myuval,mystruc,mysols] = feval(f,k,klabels,addlargs);
  while mystruc.badSSs > 5
    myc = myc+1;
    if (myc == 1) & (numel(std) > 1)
      k = kmean.*(1-std.*ones(size(kmean)));
    else
      if (numel(std) > 1)
        tmp = randn
      else
        tmp = randn(size(kmean));
      end
      k = kmean.*(1+std.*randn(size(kmean)));
    end
    [myuval,mystruc,mysols] = feval(f,k,klabels,addlargs);
  end
  mystruc.k = k;

end
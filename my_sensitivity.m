function [dudk,dutil,klabels,kmat,uvals] = my_sensitivity(k,klabels,...
  addlargs,varargin)

  %  This program can be used to estimate the sensitivity of utility
  %    values with respect to various parameter values.  After the
  %    parameter values to be tested have been determined, the program uses
  %    simple first order finite differences along each parameter direction
  %    (effectively evaluating the utility function at n+1 points where n =
  %    number of parameters to be tested) to estimate the sensitivities.
  %    It uses parallel processing so it might be a good idea to make sure
  %    the matlabworkers (see parpool/the icon at the bottom left of the
  %    "command window")
  %
  %  INPUTS
  %  k-numerical values of parameters that we wish to be tested.
  %    Sensitivity calculations will be done near this point in parameter
  %    space.
  %  klabels-the names of the parameters...necessary for many of the
  %    inflammation programs
  %  addlargs-any additional arguments that need to be handed into the
  %    utility function
  %  varargin-parameter value-name pairs used to redefine values below.
  %    parcomp: ...,'parcomp',true,... turns on parallel computation
  %    (useful for complicated models, not so much for simple ones)
  %    util: ...,'util',utility_function,... hands in a different utility
  %    function other than the parabola.  The utility function needs to
  %    take 3 inputs, k (parameter values), klabels (corresponding
  %    parameter names), and addlargs which are additional arguments/pieces
  %    of info that the function needs for its computations (e.g. measured
  %    experimental data).
  %    frac: ...,'frac',1,... will use a test spacing of 100% away (and
  %    bigger) vs the original location.  -0.5 is 50% smaller.  This allows
  %    for not-so-local sensitivity analysis (let me know if we are
  %    interested in Monte Carlo clouds around a point to randomly sample
  %    and find corresponding sensitivities near the point...it can be
  %    relatively easily done
  %
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
  
  %  Redefine default values
  for vac = 1:2:numel(varargin)
    sp.(varargin{vac}) = varargin{vac+1};
  end
  
  %%  Run the set of points in parameter space
  kmat = repmat(k,1,numel(k)+1);
  kmat(:,2:end) = kmat(:,2:end)+diag(k.*sp.frac);

  if sp.parcomp
    %  Use parallel computation
    p = gcp();
    for kc = 1:(numel(k)+1)
      my_jobs(kc) = parfeval(p,sp.util,1,kmat(:,kc),klabels,addlargs);
    end    
    uvals = zeros(1,numel(k)+1);
    for kc = 1:(numel(k)+1)
      [job_index,value] = fetchNext(my_jobs);
      uvals(job_index) = value;
      fprintf('Job #%d/%d done\n',job_index,numel(k)+1);
    end
  else
    %  Without parallel computation
    for kc = 1:(numel(k)+1)
      uvals(kc) = feval(sp.util,kmat(:,kc),klabels,addlargs);
      if isinf(uvals(kc)), keyboard; end
      fprintf('Job #%d/%d done\n',kc,numel(k)+1);
    end
  end
  
  %%  Use the results from above to calculate senstivities
  dutil = uvals(2:end)'-uvals(1);
  dudk = dutil./(k.*sp.frac);
  
  %%  Reshape everything to match the dimensions of the original k
  k = reshape(k,mk,nk);
  klabels = reshape(klabels,mk,nk);
  dudk = reshape(dudk,mk,nk);
  dutil = reshape(dutil,mk,nk);
  %  We don't worry about reshaping kmat and uvals at this time

end
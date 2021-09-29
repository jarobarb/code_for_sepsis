function [c,lsq,accept,accattempt,expected_value,nrpc,swattempt,swap] = ...
  parTemp(p,varargin)

  %  Number of runs/number of points per chain
  nrpc = 10000;
  %  Number of chains to use for parallel tempering
  Nc = 3;
  %  Number of multi-chain replicates
  Nmc = 2;
  %  Runs per swap
  rps = 10;
  %  Spit out "monitor" every simcheck runs
  simcheck = 100;
  %  Baseline epsilon/jump size
  baseline_eps = 0.009;
  %  How much we trust the measurements.  Should usually depend on standard
  %  deviations associated with the actual measurements, but for now making
  %  this value bigger trusts the measurements more and the parameter
  %  values stick closer to each other.
  mtrust = 100;
  %
  mysigj = 0.1*p;

  %local/global likelihood function exponent
  f_exp_fun = @(pd) ((9*(pd(1)-0.5)).^2+(9*(pd(2)-0.5)).^2)*...
    (((9*(pd(1)-1.5)).^2+(9*(pd(2)-1.5)).^2)+1);

  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end

  %  Indices for each chain level
  n = [1:Nc];

  %  Beta definition, can be adjusted to optimize acceptance rate
  epsilon_spread = 1./(sqrt((1/6).^(n-1)));
  epsilon = baseline_eps*epsilon_spread;
  beta = (1./epsilon_spread).^2;

  %  Store past likelihood exponentials
  fy_exp = cell(Nmc,Nc);

  % other Likelihood function: f_exp_fun = @(pd)(sum((y_meas-(pd(1)*x_meas+pd(2))).^2))/(2*sigma.^2);

  %  Actual likelihood (vs its exponent), never used
  f_exp = @(x) exp(-beta(1)*f_exp_fun(x));
  fx_exp = cell(Nmc,Nc);
  lsq = cell(Nmc,Nc);
  lsqy = cell(Nmc,Nc);
  [lsq{:}] = deal(f_exp_fun(p));
  lsqy = lsq; lsqx = lsqy;
  [fx_exp{:}] = deal(get_my_prior_exp(p,'uniform',[],[])+...
    mtrust*lsqy{1}(1));

  %initialize acceptance count to calculate acceptance rate
  M = cell(Nmc,Nc);
  [M{:}] = deal(0*p);
  accept = zeros(Nmc,Nc);
  accattempt = zeros(Nmc,Nc);

  %initialize swap count to calculate swapping rate
  swattempt = zeros(Nmc,Nc-1);
  swap = zeros(Nmc,Nc-1);

  c = cell(Nmc,Nc);
  expected_value = cell(Nmc,Nc);
  [c{:}] = deal(zeros(nrpc+1,numel(p)));
  for mcc = 1:Nmc
    for cc = 1:Nc
      c{mcc,cc}(1,:) = p;
    end
  end

  for i=1:ceil(nrpc/10)

    %parallel for loop for parallel tempering
    parfor acc=1:Nmc*Nc

      [mcc,cc] = ind2sub([Nmc,Nc],acc);

      for j = 1:rps

        %initialize id
        id = j+rps*(i-1)+1;

        %put p values into a cell and get new p values
        p = c{acc}(id-1,:);
        p_new = normrnd(p, epsilon(cc)*mysigj);

        if any(p_new <= 0)
          fy_exp{acc} = inf;
          lsqy{acc} = inf;
        else
          lsqy{acc} = f_exp_fun(p_new);
          if lsqy{acc} > 10
            lsqy{acc} = inf;
          end
          fy_exp{acc} = get_my_prior_exp(p_new,'log_uniform',[],[])+...
            mtrust*lsqy{acc};
        end

        %acceptance rule for Monte Carlo algorithm
        h = min(1,exp(-beta(cc)*(fy_exp{acc}-fx_exp{acc})));

        U = rand;
        if U < h
          %accept p_new for p and switch fx and fy likelihoods
          p = p_new;
          fx_exp{acc} = fy_exp{acc};
          accept(acc) = accept(acc)+1;
          lsqx{acc} = lsqy{acc};
        end

        %add to acceptance count
        accattempt(acc) = accattempt(acc) + 1;

        %put a,b values into cells and compute the mean
        c{acc}(id,:) = p;
        lsq{acc}(id) = lsqx{acc};
        M{acc} = M{acc}+p;
      end

    end

    %initialize id outside of the parallel for loop
    id = 10+10*(i-1)+1;

    %swapping chains to explore more space
    for mcc = 1:Nmc
      for kk = Nc:-1:2

        %from MCMC with prior code, condition for when to actually swap
        if exp((beta(kk)-beta(kk-1))*(fx_exp{mcc,kk}-fx_exp{mcc,kk-1})) > rand(1)

          %swap exponentials
          temp_exp = fx_exp{kk};
          fx_exp{kk} = fx_exp{kk-1};
          fx_exp{kk-1} = temp_exp;

          %swap variables
          temp_c = c{mcc,kk}(id,:);
          c{mcc,kk}(id,:) = c{mcc,kk-1}(id,:);
          c{mcc,kk-1}(id,:) = temp_c;

          %add to swap count
          swap(mcc,kk-1) = swap(mcc,kk-1)+1;
        end

        %try to swap chains
        swattempt(mcc,kk-1) = swattempt(mcc,kk-1)+1;

      end
    end

    if mod(i,simcheck) == 0, fprintf('%d/%d done.\n',i,ceil(nrpc/10)); end

  end

  for acc = 1:Nmc*Nc
    expected_value{acc} = M{acc}/nrpc;
  end

end

function my_prior_exp_val = get_my_prior_exp(pn,my_prior,pnorm,mysigp)

  for cc = 1:size(pn,1)
    switch my_prior
      case 'uniform'
        my_prior_exp_val(cc) = 0;
      case 'log_uniform'
        my_prior_exp_val(cc) = -sum(log(pn(cc,:)));
      case 'normal_off_center'
        my_prior_exp_val(cc) = -sum((pn(cc,:)-1.1*pnorm).^2./mysigp.^2);
      case 'normal'
        my_prior_exp_val(cc) = -sum((pn(cc,:)-pnorm).^2./mysigp.^2);
    end
  end
end
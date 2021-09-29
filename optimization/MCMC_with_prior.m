function MCMC_with_prior(varargin)

% This is markov-chain monte carlo for the lung model.  It's supposed to be
% the case that it also works for different models (by changing SGV to a
% different number) but I'm not sure whether or not that switch currently
% works.
  SGV = 11;
  global ndiv;
  %  This is the interpolant.  We initialize it to be an empty object.  It
  %  is reset to be a structure useful for interpolation later.
  global z;
  z = [];
  if SGV == 11
    ndiv = 1e10;
  end
  RD = 0;
  sim_bump = 0;
  init_param_pert = 0*0.1;
  final_fig_num = 3;
  rerun = 0;
  free_params = {'mrest'};
  save_file = 1;
  meas_window = [0,1];
  param_window_frac = [0.1,10];
  my_meas_sigm_mult = 100;
  my_para_sigm_mult = 1;
  tf = 100;
  dt = 1;
  bic = 1.2;
  cpic = 0;
  caic = 0.125;
  mic = 0;
  nic = 0;
  %number of chains
  Nc = 3;
  %  Options:  'normal', 'normal_off_center', 'uniform', 'log_uniform'
  my_prior = 'log_uniform';
  plot_every_time = 0;
  global data_amount;
  data_amount = 0;
  %number of param sets to save in one run
  N =100;
  %e.g.gibbs_randp_parallelnew(1,'Outputtest')

  %% Deal with input arguments
  if isunix
    jared_plot = 0;
  else
    jared_plot = 1;
  end
  varargin_list = '';
  noise_info = {};
  for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
    switch varargin{vac}
      case 'free_params'
        varargin_list = [varargin_list,'_free_params'];
        for fpc = 1:length(free_params)
          varargin_list = [varargin_list,'_',free_params{fpc}];
        end
      case 'noise_info'
        for nc = 1:length(noise_info)
          varargin_list = [varargin_list,'_',num2str(noise_info{nc})];
        end
      case 'meas_window'
        varargin_list = [varargin_list,'_meas_window_',...
          num2str(meas_window(1)),'_',num2str(meas_window(2))];
      case {'N','RD','sim_bump','init_param_pert','params_to_bump','final_fig_num','rerun'}

      otherwise
        varargin_list = [varargin_list,'_',varargin{vac},'_',...
          num2str(varargin{vac+1})];
    end
  end
  if ~isempty(varargin_list)
    varargin_list = ['_',varargin_list];
  end

  run_name = ['MH_MCMC_N_',num2str(N),'_rd_',num2str(RD),'_ipp_',num2str(init_param_pert),...
    varargin_list,'.mat'];

  %  Don't rerun runs that have already been run
  if (length(dir(run_name))> 0 && rerun == 0)
    final_fig_num_temp = final_fig_num;
    load(run_name);
    final_fig_num = final_fig_num_temp;
  else
    %load datafile4.mat %only need this if you don't run using batch_gibbs.m
    
    %number of parameters
    Np = length(free_params);
    %number of runs before swapping
    Nrs = 10;
    if matlabpool('size') == 0
      matlabpool 6;
    end
    
    %ranges for prior distribution [ 1/bnd to bnd ]
%     bnd = 2*ones(Np,1);
    
    %distribution of Metropolis parameters
    n = [1:Nc];
    %  Jump sizes corresponding to the variance of the transition matrix.
    %  We note that some jumps are small corresponding to local chains
    %  searching locally for minima while others are large corresponding to
    %  searching globally for minima
    if SGV == 12
      baseline_eps = 0.01*0.2;
    else
      baseline_eps = 10*0.02;
    end
    epsilon_spread = 1./(sqrt((1/3).^(n-1)));
    epsilon = baseline_eps*epsilon_spread;
    %  One over the "temperature"...higher beta means lower temperature
    %  means less jumping around, lower beta more jumping around allowed
    beta = (1./epsilon_spread).^2;
%     epsilon = (bnd/10)*0.2./sqrt(beta);
    % To increase acceptance ratio: Decrease beta or epsilon, or increase stddata.
    
    %% Obtain correct parameter values for lung:
    %  This fetches any constants we want to have for a specific model
    %  choice corresponding to SGV
    SGV_consts = get_SGV_consts(SGV,free_params);
    
    %  Get initial data and a grid (only for 2-d or 3-d probs)
    [~,mygrid] = get_initial_struc(SGV, SGV_consts, [], RD, tf);
    mygrid.real_data = RD;
    
    if isfield(mygrid,'norm_consts')
      fn = fieldnames(mygrid.norm_consts);
      for fnc = 1:length(fn)
        SGV_consts.(fn{fnc}).val = mygrid.norm_consts.(fn{fnc}).val;
        SGV_consts.(fn{fnc}).stddev = 0.1*mygrid.norm_consts.(fn{fnc}).val;
      end
      for fnc = 1:length(SGV_consts.set)
        tmp_name = SGV_consts.set(fnc).name;
        SGV_consts.set(fnc).val = SGV_consts.(tmp_name).val;
      end
      for fnc = 1:length(SGV_consts.free)
        tmp_name = SGV_consts.free(fnc).name;
        SGV_consts.free(fnc).val = SGV_consts.(tmp_name).val;
        SGV_consts.free(fnc).stddev = 0.1*...
          SGV_consts.(tmp_name).val;
      end
    end
    
    %% continue with normal code
    %baseline parameter values = A matrix values
    param_vals = SGV_consts;
    
    for fpc = 1:length(SGV_consts.free)
      free_params{fpc} = SGV_consts.free(fpc).name;
      pnorm(fpc) = SGV_consts.(free_params{fpc}).val;
    end
%     pnorm = pnorm';
    p = repmat(pnorm,Nc,1);
    noise_struc = get_noise_information(noise_info{:});
    
    %##########################################################################
    t = clock;
    Njobs = 1;
    
    %randomizing the generator
    % rand('state',sum(100*clock));
    % randn('state',sum(100*clock));
    
    accepted = ones(1,Nc);
    rejected = ones(1,Nc);
    curr_accepted = accepted;
    curr_rejected = rejected;
    adjust_jumps_num = inf;
    my_adjust_jumps_count = 0;
    desired_accratio = 0.26;
    swap = zeros(1,Nc-1);
    swattempt = 0;
    
    warning off
    
    pvec = single(zeros(Nc,Np,N));
    Evec = single(zeros(Nc,N));
    likelihood = Evec;
    prior_exp = Evec;
    
    %info = zeros(1,14);
    mysigp = my_para_sigm_mult*pnorm*noise_struc.stddev.comp_para;
    mysigm = my_meas_sigm_mult*noise_struc.stddev.data_state;
%     if SGV == 11
%       mysigm = 10000*mysigm;
%     end
    
    %  These are the standard deviations corresponding to the transition
    %  pdf/the random jumps in the mcmc
    mysigj = mysigp;
    
    param_vals.bic.val = bic;
    param_vals.cpic.val = cpic;
    param_vals.caic.val = caic;
    param_vals.mic.val = mic;
    param_vals.nic.val = nic;
    
    %randomizing initial condition    %%%%%Not sure what 'randomizing' means here
    for kk = 1:Nc
      %slow version: use for defective case
      %  Ea(kk,1) = fun_gibbs_new(p(kk,:).*pnorm);  %x is input here
      
      %exact solution version (good for all cases except defective)
      my_prior_vala_exp(kk) = get_my_prior(p(kk,:),my_prior,pnorm,mysigp);

      if SGV == 12
        [Ea(kk,1),add_info] = exact_fun_gibbs_new(p(kk,:), param_vals, ...
          data_amount, 'meas_window', meas_window, 'tspan', [0:dt:tf],...
          'stddev',mysigm);  %x is input here
      else
        if kk == 1
          Ea(kk,1) = max_likelihood_prog_w_interp(1,...
            {'mlp_arg_pairs',{'stddev',mysigm}},...
            p(kk,:), free_params);
          p0 = p;
        else
          Ea(kk,1) = Ea(1,1);
        end
      end
      pa_exp(kk) = my_prior_vala_exp(kk)-Ea(kk);
      if (kk == 1 && SGV == 12), param_vals.y_exact = add_info.y_exact; end
      
    end
    j = 0;
    
    Eb_parfor = Ea;
    pn_parfor = p;
   
    while j < N
      
      %run chain
      for ll = 1:Nrs  %run each chain Np times before swapping
        % Find a new proposal using proposal density
        pn_parfor = p.*exp(epsilon'*normrnd(0,mysigj));    %generate new set of parameters
        
        %  Keep parameter values in the parameter value window
        for pc_chain = 1:size(pn_parfor,1)
          for pc = 1:size(pn_parfor,2)
            while pn_parfor(pc_chain,pc)./p0(pc_chain,pc) < ...
                param_window_frac(1)
              pn_parfor(pc_chain,pc) = pn_parfor(pc_chain,pc)*...
                (param_window_frac(2)/param_window_frac(1));
            end
            
            while pn_parfor(pc_chain,pc)./p0(pc_chain,pc) > ...
                param_window_frac(2)
              pn_parfor(pc_chain,pc) = pn_parfor(pc_chain,pc)*...
                (param_window_frac(1)/param_window_frac(2));
            end
          end
        end

        for kk = 1:Nc                    
          %  Attach a prior to the likelihood...can be uniform or normal
          pn = pn_parfor(kk,:);
          my_prior_valb_exp = get_my_prior(pn,my_prior,pnorm,mysigp);

          if kk == 1
            for kk2 = 1:Nc
              if SGV == 12
                Eb_parfor(kk2) = exact_fun_gibbs_new(pn_parfor(kk2,:), param_vals,...
                  data_amount, 'meas_window', meas_window, 'tspan', [0:dt:tf],...
                  'stddev',mysigm);
              else
                Eb_parfor(kk2) = max_likelihood_prog_w_interp(0,...
                  {'mlp_arg_pairs',{'stddev',mysigm}},...
                  pn_parfor(kk2,:), free_params);
              end
            end
          end
          Eb = Eb_parfor(kk);
          pb_exp = my_prior_valb_exp-Eb;
          
          h = min(1,exp(beta(kk)*(pb_exp-pa_exp(kk)))*...
            (prod(pn_parfor(kk,:))./prod(p(kk,:))));
%           fprintf('%g\n',h);
          
          if h>rand(1)
            p(kk,:) = pn;
            Ea(kk) = Eb;
            pa_exp(kk) = pb_exp;
            accepted(kk) = accepted(kk) + 1;
            curr_accepted(kk) = curr_accepted(kk) + 1;
          else
            rejected(kk) = rejected(kk) + 1;
            curr_rejected(kk) = curr_rejected(kk) + 1;
          end
        end
        j = j+1;
        pvec(:,:,j) = single(p);
        Evec(:,j) = single(Ea);
        likelihood(:,j) = single(-Ea);
        prior_exp(:,j) = single(get_my_prior(p,my_prior,pnorm,mysigp));
        
        if SGV == 11
          temp_final_fig = final_fig_num;
          final_fig_num = 1000;
          if size(pvec,2) > 0
            temp_pvec_MCMC = pvec;
            pvec = pvec(:,:,1:j);
            temp_likelihood_MCMC = likelihood;
            likelihood = likelihood(:,1:j);
            temp_prior_MCMC = prior_exp;
            prior_exp = prior_exp(:,1:j);
            if plot_every_time == 1
              my_plotter;
            end
            pvec = temp_pvec_MCMC;
            likelihood = temp_likelihood_MCMC;
            prior_exp = temp_prior_MCMC;
          end
          final_fig_num = temp_final_fig;
          fprintf('%g/%g steps;  ',j,N);
          for ac = 1:length(accepted)
            fprintf('%g/',accepted(kk));
            fprintf('%g ',accepted(kk)+rejected(kk));
          end
          fprintf('accratios \n');
        end
      end
      
      %try to swap chains
      swattempt = swattempt + 1;
      for kk = Nc:-1:2
        %  kk-1 is coarser chain (larger jumps) while kk is finer chain 
        %    (smaller jumps).  If coarser chain (kk-1) has pa_exp that is
        %    closer to zero, this corresponds to a better fit wrt to data
        %    and a larger positive pa_exp(kk) value.  Hence the following
        %    inequality is always true when the coarse chain has found more
        %    likely parameter values and sometimes true when the coarse
        %    chain has worse likelihood values.  Here Eb is used as a swap
        %    variable and has no physical meaning.
        if exp((beta(kk)-beta(kk-1))*(pa_exp(kk)-pa_exp(kk-1))) > rand(1)
          Eb = Ea(kk);
          Ea(kk) = Ea(kk-1);
          pb_exp = pa_exp(kk);
          pa_exp(kk) = pa_exp(kk-1);
          Ea(kk-1) = Eb;
          pa_exp(kk-1) = pb_exp;
          pn = p(kk,:);
          p(kk,:) = p(kk-1,:);
          p(kk-1,:) = pn;
          swap(kk-1) = swap(kk-1)+1;
        end
      end
      
      if j-my_adjust_jumps_count*adjust_jumps_num > adjust_jumps_num
        my_adjust_jumps_count = floor(j/adjust_jumps_num);
        curr_accratio = curr_accepted./(curr_accepted+curr_rejected);
        fprintf('Jump sizes for each chain''s transition pdfs:  \n');
        jump_mat_str = num2str((epsilon'*mysigj)');
        for pc = 1:numel(mysigj)
          fprintf('%s:   ',jump_mat_str(pc,:));
          fprintf('%s ',free_params{pc});
          fprintf('\n');
        end
        tempsig = std(pvec(:,:,1:j),0,3);
        std_dev_str = num2str(bsxfun(@ldivide,pnorm,tempsig)');
        fprintf('Normalized standard deviation in current chain:\n');
        for pc = 1:numel(mysigj)
          fprintf('%s:   ',std_dev_str(pc,:));
          fprintf('%s ',free_params{pc});
          fprintf('\n');
          if numel(mysigj) > 1
            mysigj(pc) = max(tempsig(1,pc),(mysigp(1,pc)+tempsig(1,pc))./2);
          end
        end
        fprintf('\n');
        fprintf('Current acceptance ratio:  ');
        fprintf('%g ',curr_accratio);
        fprintf('\n');
        fprintf('%g/%g steps\n',j,N);
        baseline_eps = baseline_eps*mean(curr_accratio)/desired_accratio;
        epsilon = baseline_eps*epsilon_spread;
        curr_accepted = 0*curr_accepted;
        curr_rejected = 0*curr_rejected;
      end      
      
%         
%       if mod(j,10*length(free_params)) == 0
%         disp(sprintf('%d/%d',j,N));
%         accratio = accepted./(accepted+rejected);
%         swapratio = swap./(swap+swattempt);
%       end
      
    end
    
    accratio = accepted./(accepted+rejected);
    swapratio = swap./(swap+swattempt);
    
    %save data in a file
    save(run_name);
  end

  my_plotter;

  fprintf('accept ratios: %g;  ',accratio);
  fprintf('\n');
  fprintf('swap ratios: %g;  ',swapratio);
  fprintf('\n');
  
  if size(pvec,2) > 1
    figure(100+final_fig_num);
    temp1 = pvec(:,1,:);
    temp2 = pvec(:,2,:);
%     scatter3(temp1(:)./pnorm(1),temp2(:)./pnorm(2),exp(-Evec(:)./mysigm.^2),'*');
    plot(temp1(:)./pnorm(1),temp2(:)./pnorm(2),'.','MarkerSize',1);
    xlabel([free_params{1},' estimate/',free_params{1},' exact']);
    ylabel([free_params{2},' estimate/',free_params{2},' exact']);
  end
  
  [~,max_lik_est_ind] = max(temp_likelihood);
  
  max_lik_ps = temp1(:,:,max_lik_est_ind);
  
  max_likelihood_prog({'stddev',mysigm,'plot_comp',1},...
    max_lik_ps,free_params,'unique_identifier',...
    [2011]);
  
end

function my_prior_valb_exp = get_my_prior(pn,my_prior,pnorm,mysigp)

  for cc = 1:size(pn,1)
    switch my_prior
      case 'uniform'
        my_prior_valb_exp(cc) = 0;
      case 'log_uniform'
        my_prior_valb_exp(cc) = -sum(log(pn(cc,:)));
      case 'normal_off_center'
        my_prior_valb_exp(cc) = -sum((pn(cc,:)-1.1*pnorm).^2./mysigp.^2);
      case 'normal'
        my_prior_valb_exp(cc) = -sum((pn(cc,:)-pnorm).^2./mysigp.^2);
    end
  end
end
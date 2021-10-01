%%  %  Some initial data processing.  Need to run before making plots.
%  Running this script/pressing F5 will process all data up until the first
%  "return" below.  Data from previous figures may be needed for subsequent
%  figures, so run each figure from start to finish to see them all.
%%  This folder ignores 1.94e? data
%  The data has uncertainties associated with it making it unsuitable for
%  use in the calibrating data set
outlier_matrix_194 = [4,1,1;4,1,2;4,1,3;4,1,4;4,1,5;4,1,6;...
  4,2,1;4,2,2;4,2,3;4,2,4; 4,3,1;4,3,2;4,4,1];
tmptd = 24;
plot_all_data(outlier_matrix_194,'td',tmptd)

%%  Location of paths (as of 2019/05/09)
addpath ./optimization/
addpath ./Finalized_Code/
%  This is for comparing structures...otherwise, not needed.
addpath ./Matlab_tools/

%%  List of optimal parameter sets
%  Use the last one

%%  What we need for a successful sensitivity analysis-V1
klabelsdef = {'B_infy','f','kD1','Bmf','nu4'};
kcurr = [5.0046602934939192764
         70.999380840544219495
       0.030009711419666851295
         80.838762887398843304
         1468.9689068779709942]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0],'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2'}};
%  Defines our "cost function"
fdebug_script;

%%  Original values from V1 push td back to 24
klabelsdef = {'B_infy','f','kD1','Bmf','nu4'};
kcurr = [4.98880494497202
          65.8337684515266
        0.0151293905016911
          46.3653249934213
          2587.35329191681]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,24],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','td'}};
fdebug_script;

%%  Original values from V1+k1 push td back to 24
klabelsdef = {'B_infy','f','kD1','Bmf','nu4','k1'};
kcurr = [0.681953653737825
          199.947235257306
        0.0141469258709592
          820.033777707272
          4001.52197195351
           1.3127433029186]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,600,0,0,24],...
  'fixedklabels',{'use_log_Bc','clot_capacity','kD3','kD2','td'}};
fdebug_script;

%%  V1 plus clot_capacity
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4'};
kcurr = [4.99749753688595
          68.5727998150189
        0.0236375115255366
          66.9862191611089
           379.28607637041
          1569.42348460266]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],'fixedklabels',{'use_log_Bc','kD3','kD2'}};
%  Defines our "cost function"
fdebug_script;

%%  V1 plus clot capacity, k1, and mul, minus nu4
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','mul','k1'};
kcurr = [0.728530455732297
          214.397633476279
        0.0327149818566835
          1217.57723740552
          422.227256349954
       0.00200167871174815
          1.27747099703467]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,3640.20598805358],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','nu4'}};
fdebug_script;

%%
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','mul','k1'};
kcurr = [0.728530455732297
          214.397633476279
        0.0327149818566835
          1217.57723740552
          422.227256349954
       0.00200167871174815
          1.27747099703467]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',false,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0,3640.20598805358],...
  'fixedklabels',{'use_log_Bc','kD3','kD2','nu4'}};
usegenalg = false; fdebug_script;

%%  V1+clot_capacity+k1
klabelsdef = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
kcurr = [0.736964699290478
          213.335614287955
        0.0344697636772473
          1207.32981781487
          407.214725448303
          3640.09871963404
          1.27176887084432]';
addlargs = {'penalty_vector',[0 0 0 0],'log_transformed',true,...
  'categorical_penalty',{{128,0},{248,1},{505,2}},'use_all',true,...
  'solve_for_kA',true,'outlier_matrix',outlier_matrix_194,'fixedk',...
  [0,0,0],'fixedklabels',{'use_log_Bc','kD3','kD2'}};
fdebug_script;

%%  Find/deal with kA and fixed (run only once per above parameter values)
%  We pretty much always solve for kA so that way anti-inflammatory
%  inhibition is 75%:  1/(1+kA*Amax) = 0.75.  Thus we rarely know what kA
%  is as it is actually changing to make sure 1/(1+kA*Amax) is always 0.75
%  = inhi_perc (which theoretically we can change if we want, but we
%  haven't explored other values of inhi_perc at this time).
[my_quantity,stats,sol] = fdebugvarg(kcurr,'plot_sol',false);
klabelsdef = {klabelsdef{:},'kA'};
kcurr(end+1) = sol{3}.struc.fp.kA;
for ac = 1:numel(addlargs)
  if isequal(addlargs{ac},'solve_for_kA')
    addlargs{ac+1} = false;
  end
end

%  Get names of all parameters (and their descriptions)
for aac = 1:numel(addlargs)
  if isequal(addlargs{aac},'fixedk')
    kcurr = [kcurr,addlargs{aac+1}];
  elseif isequal(addlargs{aac},'fixedklabels')
    klabelsdef = {klabelsdef{:},addlargs{aac+1}{:}}
  end
end

fdebug_script;

%%  Sensitivity estimates involving all parameters (even ones not in use)
%  Retrieves utility and measurement based relative sensitivities for all
%  parameters including ones that aren't in use or ones that are logicals
%  (e.g. true/false).  Many of these we aren't interested in, but this
%  initial sensitivity analysis helps us get rid of those.  Because many of
%  the parameters aren't in use, the sensitivity matrix is singular.
[~,~,~,~,kall,klabelsall,desc] = get_parameters(kcurr,klabelsdef);

%  Run the sensitivity analysis
[dudkall,dutilall,klabelsoutall,kmatall,uvalsall,s_strucall] = ...
  my_sensitivity_nongeneral(kall,klabelsall,addlargs,'util',reg_util,...
  'parcomp',true,'frac',sqrt(eps),'std',0);

%%  Loop through all variables individually
%  Different way to find non-used parameters/logicals
for kc = 1:numel(kall)
  inds = ismember(klabelsall,klabelsall{kc});
  [dudk,dutil,klabelsout,kmat,uvals,s_struc] = ...
    my_sensitivity_nongeneral(kall(inds),klabelsall(inds),...
    {addlargs{:},'fixedk',kall(~inds),'fixedklabels',{klabelsall{~inds}}},...
    'util',reg_util,'parcomp',true,'frac',0.01,'npts',18);
  fprintf('%10g %s\n',dudk,klabelsall{kc});
end

%%  Print out the sensitivities (many which will be 0/nan)
[~,inds] = sort(abs(s_strucall.smsqr));
for j = 1:numel(kall)
  ind = inds(j);
  name = klabelsall{ind};
  if isfield(desc.alt_name,name), name = desc.alt_name.(name); end
  fprintf('%s,%3.6e,%g,%g,%s\n',name,kall(ind),...
    s_strucall.smsqr(ind),dudkall(ind),desc.descs.(klabelsall{ind}));
end

%%  Take out parameters that don't seem to matter
inds = ((abs(s_strucall.smsqr') > 1e-10) &...
  (~isnan(s_strucall.smsqr'))) | ismember(klabelsall,{'nu4'});
% inds = inds & ~eps0ind & ~gammaind;

klabelsimp = {klabelsall{inds}};
kimp = kall(inds);
for kc = 1:numel(kcurr)
  tmp = ismember(klabelsimp,klabelsdef{kc});
  if ~any(tmp)
    klabelsimp = {klabelsimp{:},klabelsdef{kc}};
    kimp(end+1) = kcurr(kc);
  end
end

%  Take out any unwanted guys:
inds = ismember(klabelsimp,{'tscale','c1','use_log_Bc','kD3','kD2','td'});
klabelsimp = {klabelsimp{~inds}};
kimp = kimp(~inds);


%%  Sensitivities only when desired behavior is produced
%  That is sensitivity only for used non-logical parameters
[dudk,dutil,klabelsout,kmat,uvals,s_struc] = ...
  my_sensitivity_nongeneral(kimp,klabelsimp,addlargs,'util',reg_util,...
  'parcomp',true,'frac',sqrt(eps),'std',0*sqrt(eps));

%%  Test just 1 (or a handful) of the variables
% Skip this if you want all results
inds = ismember(klabelsimp,{'kBl','sl'});
% inds = ismember(klabelsimp,{'kBl','sl'});
[dudkind,dutilind,klabelsoutind,kmatind,uvalsind,s_strucind] = ...
  my_sensitivity_nongeneral(kimp(inds),klabelsimp(inds),...
  {addlargs{:},'fixedk',kimp(~inds),'fixedklabels',{klabelsimp{~inds}}},...
  'util',reg_util,'parcomp',false,'frac',sqrt(eps),'std',0);

%%  Relative sensitivities (using all "important" variables)
fprintf('Sensitivities list:\n');
[~,inds] = sort(abs(s_struc.smsqr));
% inds = 1:numel(kimp);
for j = 1:numel(kimp)
  ind = inds(j);
  name = klabelsimp{ind};
  if isfield(desc.alt_name,name), name = desc.alt_name.(name); end
  fprintf('%s,%3.3g,%3.3g,%g,%s\n',name,kimp(ind),...
    s_struc.smsqr(ind),dudk(ind),desc.descs.(klabelsimp{ind}));
end

%%  Same list reformatted
inds = 1:numel(kimp);
for j = 1:numel(kimp)
  ind = inds(j);
  name = klabelsimp{ind};
  if isfield(desc.alt_name,name), name = desc.alt_name.(name); end
  fprintf('%s,%3.3g\n',name,kimp(ind));
end

%%  Plot for some of the sensitivities...not that helpful
%  Theoretically, rows that are the same are possibly collinear.  In this
%  plot, however, a select few squares dominate making it hard to compare
%  rows.
close(figure(373)); figure(373)
pcolor(log(abs(s_struc.dsmmat))); colorbar

%%  Store shifted simulated measurements and parameter values
%  m parameters by n observations
dsmmat = s_struc.dsmmat;
smmat_shif = s_struc.smmat_shif;
kmat_shif = s_struc.kmat_shif(2:end,:);

%%  Calculate collinearities/fit quality for a bunch
maxcombos = 7;
clear cl;

for cc = 2:maxcombos
  cl(cc).C = nchoosek([1:size(dsmmat,1)],cc);
  for Cc = 1:size(cl(cc).C,1)
    Sc = dsmmat(cl(cc).C(Cc,:),:);
    lambdas = svd(Sc*Sc');
    cl(cc).CI(Cc) = 1./sqrt(min(lambdas));
    
    coefs = (kmat_shif(cl(cc).C(Cc,:),:)')\(smmat_shif);
    SSE = sum(sum(((kmat_shif(cl(cc).C(Cc,:),:)')*coefs-smmat_shif).^2));
    SST = sum(sum(smmat_shif.^2));
    cl(cc).R2(Cc) = 1-SSE/SST;
    cl(cc).ss(Cc) = sum(s_struc.smsqr(cl(cc).C(Cc,:)));
    if mod(Cc,10000) == 0
      fprintf('%d: %g/%g\n',cc,Cc,size(cl(cc).C,1)); 
    end
  end
end

%%  Display best candidates
nbest = 10;
%  Display n best candidates info based on R2 values
for cc = 2:maxcombos
  [a,b] = sort(cl(cc).R2);
  fprintf('%d:\n',cc);
  for j = nbest:-1:1
    ind = b(end-(j-1));
    fprintf('%10g',cl(cc).R2(ind));
    fprintf('%10g',cl(cc).ss(ind));
    fprintf('%10g',cl(cc).CI(ind));
    fprintf('%10s',klabelsimp{cl(cc).C(ind,:)});
    fprintf('\n');
  end
end

%  Display n best candidates info based on CI values
for cc = 2:maxcombos
  [a,b] = sort(cl(cc).CI,'descend');
  fprintf('%d:\n',cc);
  for j = nbest:-1:1
    ind = b(end-(j-1));
    fprintf('%10g',cl(cc).CI(ind));
    fprintf('%10g',cl(cc).ss(ind));
    fprintf('%10g',cl(cc).R2(ind));
    fprintf('%10s',klabelsimp{cl(cc).C(ind,:)});
    fprintf('\n');
  end
end

%  Display n best candidates info based on CI values
for cc = 2:maxcombos
  [a,b] = sort(cl(cc).ss,'descend');
  fprintf('%d:\n',cc);
  for j = nbest:-1:1
    ind = b(end-(j-1));
    fprintf('%10g',cl(cc).ss(ind));
    fprintf('%10g',cl(cc).R2(ind));
    fprintf('%10g',cl(cc).CI(ind));
    fprintf('%10s',klabelsimp{cl(cc).C(ind,:)});
    fprintf('\n');
  end
end

%  Displan n best candidates based on R2 values with a CI value < CI_crit
CI_crit = 25;
for cc = 2:maxcombos
  [a,b] = sort(cl(cc).R2.*(cl(cc).CI<CI_crit));
  fprintf('%d:\n',cc);
  for j = nbest:-1:1
    ind = b(end-(j-1));
    fprintf('%10g',cl(cc).R2(ind));
    fprintf('%10g',cl(cc).ss(ind));
    fprintf('%10g',cl(cc).CI(ind));
    fprintf('%10s',klabelsimp{cl(cc).C(ind,:)});
    fprintf('\n');
  end
end

return

%%  Calculate one collinearity index/fitting capability
inds = find(ismember(klabelsimp,klabelsdef));
% inds = find(ismember(klabelsimp,...
%   {'kD1','Bmf','clot_capacity','k1','kBl','mul','sA'}));

Sc = dsmmat(inds,:);
lambdas = svd(Sc*Sc');
ColInd = 1./sqrt(min(lambdas))

coefs = (kmat_shif(inds,:)')\smmat_shif;
SSE = sum(sum(((kmat_shif(inds,:)')*coefs-smmat_shif).^2));
SST = sum(sum(smmat_shif.^2));
R2 = 1-SSE/SST

ss = sum(s_struc.smsqr(inds))

%  This local fit doesn't work because kBl and sl are collinear (see below)
%  hence not locally identifiable
[beta,Sigma,E,CovB,logL] = mvregress(kmat_shif(inds,:)',smmat_shif);

%%  Calculate one collinearity index/fitting capability
inds = find(ismember(klabelsimp,{'kBl','sl'}));

Sc = dsmmat(inds,:);
lambdas = svd(Sc*Sc');
ColInd = 1./sqrt(min(lambdas))

coefs = (kmat_shif(inds,:)')\smmat_shif;
SSE = sum(sum(((kmat_shif(inds,:)')*coefs-smmat_shif).^2));
SST = sum(sum(smmat_shif.^2));
R2 = 1-SSE/SST

[beta,Sigma,E,CovB,logL] = mvregress(kmat_shif(inds,:)',smmat_shif);
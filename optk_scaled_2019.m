function [my_quantity,stats,sol] = optk_scaled_2019(k,klabels,varargin)

%all concentrations are in 10^6 units
plot_solution = false;
penalty_vector = [0 0 0 0];
%  Note that regularization takes place in logarithmic space.  See below.
regularization_coefs = 0*k;
regularization_lower = zeros(size(k));
regularization_upper = inf(size(k));
categorical_penalty = {};
%  Rescales the measured bacterial data if set to true (alternatively, one
%  can think of this as rescaling the simulated bacterial data since
%  SSE/SST is unaffected if one chooses to measure SSE and SST in "measured
%  units" vs "simulated units"
log_transformed = false;
use_all = false;
plot_sol = false;
parcomp = true;
eval_stats = false;
fig_bump = 0;
save_for_debug = false;
solve_for_kA = false;
fixedk = [];
fixedklabels = {};
plot_b_terms = true;

mycols = 'rbgk';
%  Picking out outliers is involved...
%    outlier_matrix(i,1) is the index corresponding to the bacterial load
%    for a give observation (128 => outlier(i,1) = 1, 248 => outlier(i,1) =
%    2, etc).  outlier_matrix(i,2) is the index corresponding to the
%    experiment's index for that particular bacterial load.
%    outlier_matrix(i,3) is the index corresponding to the time at which
%    the measurement was made
outlier_matrix = [];
classify_vec = false(size(penalty_vector));

%  varargin stuff
for vac = 1:2:numel(varargin)
  eval([varargin{vac},' = varargin{vac+1};']);
end

if (numel(regularization_coefs) ~= numel(k)) | ...
    (numel(regularization_lower) ~= numel(k)) | ...
    (numel(regularization_upper) ~= numel(k))
  error('regularization_stuff is the wrong size');
else 
  %  Makes sure regularization stuff has right dimensions
  tmp = regularization_coefs;
  regularization_coefs = k;
  regularization_coefs(:) = tmp(:);
  tmp = regularization_lower;
  regularization_lower = zeros(size(k));
  regularization_lower(:) = tmp(:);
  tmp = regularization_upper;
  regularization_upper = zeros(size(k));
  regularization_upper(:) = tmp(:);
end

%  Must extract "td" if it exists
tmptd = 0;
for vac = 1:2:numel(varargin)
  if isequal(varargin{vac},'td')
    tmptd = varargin{vac+1};
  elseif isequal(varargin{vac},'fixedk')
    tmpfixedk = varargin{vac+1};
  elseif isequal(varargin{vac},'fixedklabels')
    tmpfixedklabels = varargin{vac+1};
  end
end
for kc = 1:numel(klabels)
  if isequal(klabels{kc},'td')
    tmptd = k(kc);
  end
end
if exist('tmpfixedk')
  for kc = 1:numel(tmpfixedklabels)
    if isequal(tmpfixedklabels{kc},'td')
      tmptd = tmpfixedk(kc);
    end
  end
end
[exd] = get_data(outlier_matrix,'td',tmptd);

% fprintf('td = %g\n',min(exd.t.list));

% %  Added back in 8h info
% time4 = [2 4 5 6 8 48]; %Corresponds with Bsource4
% actual_conc4 = [312.5 150 50 190 800 145]; %averages of the concentrations %removed time=8, 800
% 
% %  Added back in 4h info
% time3 = [0 2 2.5 3 4 6 8 18 28]; %Corresponds with Bsource3
% actual_conc3 = [170 170 144 50 495 150 190 50 50]; %averages of the concentrations %removed time=4, 495
% 
% time3 = [2 2.5 3 6 8 18 28]; %Corresponds with Bsource3
% actual_conc3 = [170 144 50 150 190 50 50]; %averages of the concentrations %removed time=4, 495
% 
% time2 = [4 8 24 48 72]; %Corresponds with Bsource2
% actual_conc2 = [150 144 190 8 2.4]; %averages of the concentrations
% 
% time1 = [24 48 72 96]; %Corresponds with Bsource1
% actual_conc1 = [35.35 52.15 1.97 .685]; %averages of the concentrations

%plot(time1,actual_conc1,'b')
%ylabel('bacteria in blood')
%xlabel('time')
%hold on
%plot(time2,actual_conc2,'g')
%plot(time3,actual_conc3,'r')
%plot(time4,actual_conc4,'k')
%legend('Bsource=1.28e8','Bsource=2.48e8','Bsource=5.05e8','Bsource=1.94e9')

%original concentrations
% Bsource1 = 128;
% Bsource2 = 248;
% Bsource3 = 505;
% Bsource4 = 1940;
Bsourcee1 = 129;
Bsourcee2 = 247;
Bsourceem = (exd.Bsource(1)+exd.Bsource(2))/2;

bottom_guess = 0;
top_guess = max(exd.Bsource);

%  If we are planning on possibly imposing a categorical penalty, we adjust
%  the "classify_vec" so that way model_scaled_2018 gets the classification
%  (septic, aseptic, healthy) in addition to normal integration/simulation.
for cpc = 1:numel(categorical_penalty)
  if numel(categorical_penalty{cpc}{1}) == 1
    if any(categorical_penalty{cpc}{1} == exd.Bsource)
      ind = find(categorical_penalty{cpc}{1} == exd.Bsource);
      classify_vec(ind) = true;
    end
  end
end

[~,sim.B.ave,sol,rp,op,fp,strp,kcellfull,klabelcellfull] = model_scaled_2018_loop([],[],...
  exd,k,klabels,'plot_solution',plot_solution,'parcomp',parcomp,...
  'outlier_matrix',outlier_matrix,'classify_vec',classify_vec,...
  'solve_for_kA',solve_for_kA,'fixedk',fixedk,'fixedklabels',fixedklabels,...
  'plot_b_terms',plot_b_terms);
%  kA check---informal estimate.  Only exact if worst case is septic and we
%  had a constraint that required its classification as such.
if strp{1}.check_kA
  for sc = 1:numel(sol)
    tmpA(sc) = sol{sc}.anti_inf(end);
  end
  fprintf('1/(Asep*kA+1) = approx %g\n',1/(max(tmpA)*fp{3}.kA+1));
end

%  Finds simulation values corresponding to all measurements (except
%  outliers) instead of just 1 simulated value per "average" measurement.
sim.B.list = [];
%  This needs to copy get_data exactly
for Bloadc = 1:numel(exd.t.all)
  for expc = 1:numel(exd.t.all{Bloadc})
    for tc = 1:numel(exd.t.all{Bloadc}{expc})
      t = exd.t.all{Bloadc}{expc}(tc);
      sim.B.all{Bloadc}{expc}(tc) = interp1(exd.t.ave{Bloadc},...
        sim.B.ave{Bloadc},t,'linear','extrap');
      if ~exd.outlier{Bloadc}{expc}(tc)
        sim.B.list(end+1) = sim.B.all{Bloadc}{expc}(tc);
      end
    end
  end
end

%Use this method
Bstave = [[exd.Bs.ave{:}]',[exd.t.ave{:}]'];
Bave = [exd.B.ave{:}]';
Bsimave = [sim.B.ave{:}]';
logindsall = sim.B.list > 0;
logindsave = Bsimave > 0;
logBsimall = log(sim.B.list(logindsall));
logBlist = log(exd.B.list(logindsall));
logBsimave = log(Bsimave(logindsave));
logBave = log(Bave(logindsave));
% 
% Bmfall = sum(exd.B.list.*sim.B.list)./sum(sim.B.list.^2);
% Bmfave = sum(Bave.*Bsimave)./sum(Bsimave.^2);
% Bmflogall = exp(sum(logBlist-logBsimall)./numel(logBsimall));
% Bmflogave = exp(sum(logBave-logBsimave)./numel(logBsimave));
% 
% if nonunity_Bmf
%   if use_all
%     if log_transformed
%       Bmf = Bmflogall;
%     else
%       Bmf = Bmfall;
%     end
%   else
%     if log_transformed
%       Bmf = Bmflogave;
%     else
%       Bmf = Bmfave;
%     end
%   end
% else
%   Bmf = 1; Bmflogall = 1; Bmfall = 1; Bmflogave = 1; Bmfave = 1;
% end

stats.simallmeas = sim.B.list;
stats.simavemeas = Bsimave;
stats.simlogallmeas = logBsimall;
stats.simlogavemeas = logBsimave;
stats.expallmeas = exd.B.list;
stats.expavemeas = Bave;
stats.explogallmeas = logBlist;
stats.explogavemeas = logBave;

stats.log_transformed = log_transformed;
stats.use_all = use_all;
if log_transformed
  if use_all
    stats.simmeas = stats.simlogallmeas; 
    stats.expmeas = stats.explogallmeas;
  else
    stats.simmeas = stats.simlogavemeas; 
    stats.expmeas = stats.explogavemeas;
  end
else
  if use_all
    stats.simmeas = stats.simallmeas; 
    stats.expmeas = stats.expallmeas;
  else
    stats.simmeas = stats.simavemeas; 
    stats.expmeas = stats.expavemeas;
  end
end

stats.SSEall = sum((exd.B.list-sim.B.list).^2);
stats.SSTall = sum((exd.B.list-mean(exd.B.list)).^2);
stats.SSEave = sum((Bave-Bsimave).^2);
stats.SSTave = sum((Bave-mean(Bave)).^2);
stats.SSElogall = sum((logBlist-logBsimall).^2);
stats.SSTlogall = sum((logBlist-mean(logBlist)).^2);
stats.SSElogave = sum((logBave-logBsimave).^2);
stats.SSTlogave = sum((logBave-mean(logBave)).^2);

% if isinparallel | (~parcomp)
%   for Bloadc = 1:numel(exd.Bsource)
%     [sepsis_exp{Bloadc},sol{Bloadc}] = model_scaled_2018('k',k,...
%       'klabels',klabels,'time',exd.t.ave{Bloadc},...
%       'Bsource_in',exd.Bsource(Bloadc),...
%       'plot_solution',plot_solution,'parcomp',parcomp);
%   end
% else
%   parfor Bloadc = 1:numel(exd.Bsource)
%     [sepsis_exp{Bloadc},sol{Bloadc}] = model_scaled_2018('k',k,...
%       'klabels',klabels,'time',exd.t.ave{Bloadc},...
%       'Bsource_in',exd.Bsource(Bloadc),...
%       'plot_solution',plot_solution,'parcomp',parcomp);
%   end  
% end

val2 = 0;
if any(penalty_vector)
  kcell = num2cell(k);
  mycell = {klabels{:};kcell{:}};
  vals = my_utility2(mycell{:});
  %  A penalty to try and make the bacterial levels for the healthy
  %  situation nonzero (cause they usually are super close to zero)
  val2 = penalty_vector(1)*vals(1)+penalty_vector(2)*vals(2)+...
    penalty_vector(3)*vals(3)+...
    penalty_vector(4)*(max(sol{1}.bacteria)-max(exd.B.ave{1}))^2/...
    max(exd.B.ave{1})^2;
%   val2 = val2+(sepsis_exp1(2)-exd.B.ave{1}(2))^2;
  
  %  These guys return 
%   f1 = my_finder(mycell{:},'Bsource_in',Bsource1);
%   f2 = my_finder(mycell{:},'Bsource_in',Bsource2);
  f1 = 0; f2 = 2;
  % fe1 = my_finder(mycell{:},'Bsource_in',Bsourcee1);
  % fe2 = my_finder(mycell{:},'Bsource_in',Bsourcee2);
  % fem = my_finder(mycell{:},'Bsource_in',Bsourceem);
  % [sepsis_expe1,sole1]=model_scaled_2018(k,klabels,[],[],...
  %   Bsourcee1,'plot_solution',plot_solution);
  % [sepsis_expe2,sole2]=model_scaled_2018(k,klabels,[],[],...
  %   Bsourcee2,'plot_solution',plot_solution);
else
  f1 = 0; f2 = 0;
end
for rc = 1:numel(k)
  if k(rc) < regularization_lower(rc)
    val2 = val2+regularization_coefs(rc).*...
      (log(regularization_lower(rc))-log(max(abs(k(rc)),eps(abs(k(rc)))))).^2;
  end
  if k(rc) > regularization_upper(rc)
    val2 = val2+regularization_coefs(rc).*...
      (log(regularization_upper(rc))-log(max(abs(k(rc)),eps(abs(k(rc)))))).^2;
  end
end

%  Code added 2020/09/25.  Tries to make sure septic simulations are still
%  high and aseptic simulations still low by the end of the simulation.
%  Note that for fint = inf, this is not "infty" but rather a really large
%  number depending on how long our integrator integrates "until steady
%  state".
for myc = 1:numel(sol)
  if isinf(fp{1}.fint)
    yfin(myc) = sol{myc}.bacteria(end);
  else
    yfin(myc) = interp1(sol{myc}.solution_time,sol{myc}.bacteria,...
      fp{1}.fint,'linear','extrap');
  end
end

nondesired_outcomes = 0;
outcomex = [];
outcomey = [];

if numel(categorical_penalty) > 0
  kcell = num2cell(k);
  mycell = {klabels{:};kcell{:}};
  for cpc = 1:numel(categorical_penalty)
    if numel(categorical_penalty{cpc}{1}) == 1
      outcomex(end+1) = categorical_penalty{cpc}{1};
      if any(categorical_penalty{cpc}{1} == exd.Bsource)
        ind = find(categorical_penalty{cpc}{1} == exd.Bsource);
        outcomey(end+1) = sol{ind}.classification;
      else
        outcomey(end+1) = my_finder(mycell{:},'Bsource_in',outcomex(end));
      end
      if ~ismember(outcomey(end),categorical_penalty{cpc}{2})
        nondesired_outcomes = nondesired_outcomes+10;
      end
      if ismember(1,categorical_penalty{cpc}{2})
%         if fp{cpc}.Bfinasepmax < 1/eps
          nondesired_outcomes = nondesired_outcomes+...
            0.1*(yfin(cpc) > fp{cpc}.Bmax*fp{cpc}.Bfinasepmax).*...
            (yfin(cpc) - fp{cpc}.Bmax*fp{cpc}.Bfinasepmax).^2;
%         end
      elseif ismember(2,categorical_penalty{cpc}{2})
%         if fp{cpc}.Bfinsepmin > eps
          nondesired_outcomes = nondesired_outcomes+...
            0.1*(yfin(cpc) < fp{cpc}.Bmax*fp{cpc}.Bfinsepmin).*...
            (yfin(cpc) - fp{cpc}.Bmax*fp{cpc}.Bfinsepmin).^2;
%         end
      end
    else
      nondesired_outcomes = nondesired_outcomes+10;
      lowx = categorical_penalty{cpc}{1}(1);
      higx = categorical_penalty{cpc}{1}(2);
      tolx = categorical_penalty{cpc}{1}(3);
      desval = categorical_penalty{cpc}{2};
      kx = find(lowx == outcomex);
      if isempty(kx)
        lowval = my_finder(mycell{:},'Bsource_in',lowx);
      else
        lowval = outcomey(kx);
      end
      kx = find(higx == outcomex);
      if isempty(kx)
        higval = my_finder(mycell{:},'Bsource_in',higx);
      else
        higval = outcomey(kx);
      end
      if (lowval == tolx) | (higval == tolx)
        nondesired_outcomes = nondesired_outcomes-10;
      %  Assume higval > lowval...i.e. "f" is monotonically increasing
      elseif (higval < desval) | (lowval > desval)
        %  Do nothing
      else
        %  Very specific to values being either 0, 1, or 2
        while higx-lowx > tolx
          midx = (lowx+higx)/2;
          midval = my_finder(mycell{:},'Bsource_in',midx);
          if midval == categorical_penalty{cpc}{2}
            nondesired_outcomes = nondesired_outcomes-10;
            break;
          elseif midval == lowval
            lowx = midx;
          else
            higx = midx;
          end
        end
      end
    end
  end
%   fprintf('\n');
end

%{
% calculated_conc = bacteria_blood;
%we aren't using this method
my_quantity1 = (sum((actual_conc1 - sepsis_exp1).^2))/(sum(actual_conc1.^2));
my_quantity2 = (sum((actual_conc2 - sepsis_exp2).^2))/(sum(actual_conc2.^2));
my_quantity3 = (sum((actual_conc3 - sepsis_exp3).^2))/(sum(actual_conc3.^2));
my_quantity4 = (sum((actual_conc4 - sepsis_exp4).^2))/(sum(actual_conc4.^2));

my_quantity_old = my_quantity1 +  my_quantity2 +  my_quantity3 + my_quantity4;
%}
%  The more natural/easier to explain least squares quantity, sum of
%  squared errors (SSE) for all pieces of data predicted divided by the sum
%  of squared totals (SST) for the actual data.  Note that R^2 (coefficient
%  of multiple determination) = 1-SSE/SST and that R^2 gives the amount of
%  variance "explained" by the model.

% Bmfall1 = fminsearch(@(Bmf) sum((exd.B.list-Bmf*sim.B.list).^2),1)
% Bmfave1 = fminsearch(@(Bmf) sum((Bave-Bmf*Bsimave).^2),1)
% Bmflogall1 = fminsearch(@(Bmf) ...
%   sum((logBlist-log(Bmf*exp(logsimall))).^2),1)
% Bmflogave1 = fminsearch(@(Bmf) ...
%   sum((logBave-log(Bmf*exp(logBsimave))).^2),1)

%  Too many options...not currently included
% for Bloadc = 1:numel(exd.Bsource)
%   SSEc{Bloadc} = sum((exd.B.ave{Bloadc}-Bmfave*sepsis_exp{Bloadc}).^2);
%   SSTc{Bloadc} = sum((exd.B.ave{Bloadc}-...
%     mean([exd.B.ave{:}])).^2);
%   SSElc{Bloadc} = sum((log(exd.B.ave{Bloadc})-...
%     Bmflogave*log(sepsis_exp{Bloadc})).^2);
%   SSTlc{Bloadc} = sum((log(exd.B.ave{Bloadc})-...
%     mean(log(Bave(inds)))).^2);
%   SSEallc{Bloadc} = sum((exd.B.list-sim.B.list).^2);
%   SSTallc{Bloadc} = sum((exd.B.list-mean(exd.B.list)).^2);
% end

% SSEall = {}; SSTall = {};
% for Bloadc = 1:numel(time)
%   for expc = 1:numel(exd.t.all{Bloadc})
%     for tc = 1:numel(exd.t.all{Bloadc}{expc})
%       if ~exd.outlier{Bloadc}{expc}(tc)
%         ind = find(exd.t.all{Bloadc}{expc}(tc) == tunique);
%         actual_concs{Bloadc}(ind) = actual_concs{Bloadc}(ind)+...
%           exd.B.all{Bloadc}{expc}(tc);
%         ns{Bloadc}(ind) = ns{Bloadc}(ind)+1;
%       end
%     end
%   end
%   exd.B.ave{Bloadc} = actual_concs{Bloadc}./ns{Bloadc};
% end
if use_all
  stats.SSEl = stats.SSElogall;
  stats.SSTl = stats.SSTlogall;
  stats.SSE = stats.SSEall;
  stats.SST = stats.SSTall;
else
  stats.SSEl = stats.SSElogave;
  stats.SSTl = stats.SSTlogave;
  stats.SSE = stats.SSEave;
  stats.SST = stats.SSTave;
end

%SSE=SSE2;
%SST=SST2;
%  Penalties try to make it so that the maximum Bsource corresponds to
%  septic death, the minimum Bsource corresponds to a healthy outcome,
%  and at least one of the middle guys is an aseptic death.
% my_quantity = (SSE/SST)+1000*include_penalties*...
%   (~(sol4.bacteria(end)>1)+~(sol1.damage(end)<1)+...
%   exp(-max(sol1.bacteria)/100)+...
%   0*(~(((sol2.bacteria(end)<1) & (sol2.pro_inf(end) > 0.5)) |...
%   ((sol3.bacteria(end)<1) & (sol3.pro_inf(end) > 0.5)))));
% my_quantity = 1e6*(sol1.damage(end)<1)-max(sol1.bacteria);
% my_quantity = (SSE/SST)+1000*include_penalties*...
%   (~(sol4.bacteria(end)>1)+~(sol3.bacteria(end)>1)+...
%   ~(sol3.bacteria(end)>1)+~(sol1.damage(end)<1)+...
%   ~((sole1.bacteria(end)<1)&(sole1.damage(end)>1))+...
%   ~((sole2.bacteria(end)<1)&(sole2.damage(end)>1)));

%  If "clot_capacity" is larger than 1940, it has no effect on results and
%  is thus not uniquely determined.  To uniquely determine it, we add a
%  one-sided quadratic dependence (if it is among the set of parameters to
%  be fit).
ind = find(ismember(klabels,'clot_capacity'));
if ~isempty(ind)
  val2 = val2+(k(ind)-1940)^2*(k(ind)>1940);
end
%  If [fM-T] is always < 0, then the corresponding term in the damage
%  equation has no effect on results and in those cases f and T are not
%  uniquely determined.  To uniquely determine it, we add a one-sided
%  quadratic dependence (if it is among the set of parameters to be fit).
indf = find(ismember(klabels,'f'));
indT = find(ismember(klabels,'T'));
if isempty(indf)
  tmpf = fp{1}.f;
else
  tmpf = k(indf);
end
if isempty(indT)
  tmpT = op{1}.T;
else
  tmpT = k(indT);
end
if (~isempty(indf)) | (~isempty(indT))
  maxM = 0;
  for mMc = 1:numel(sol)
    maxM = max([maxM,max(sol{mMc}.pro_inf)]);
  end
  tmpquant = tmpf*maxM-tmpT;
  val2 = val2+tmpquant^2.*(tmpquant < 0);
end

%  Various quantities to maximize...original
% my_quantity = (SSE/SST)+1000*penalty_vector*...
%   (~(f1 == 0)+~(f2 == 2))+val2;
%  No data, just qualitative agreement with desired behavior:
stats.regularizations = val2;
stats.badSSs = nondesired_outcomes;
stats.curragreement = log_transformed*(stats.SSEl/stats.SSTl)+...
  (~log_transformed)*(stats.SSE/stats.SST);
my_quantity = log_transformed*(stats.SSEl/stats.SSTl)+...
  (~log_transformed)*(stats.SSE/stats.SST)+val2+nondesired_outcomes;
stats.my_quantity = my_quantity;

% my_quantity = SSE
% my_quantity = (SSE/SST)+val2;
%  Use of varargin allows us to plot the data for comparison

  if plot_sol
    close all
    kcell = num2cell(k);
%     klabels{end+1} = 'solve_for_kA'; kcell{end+1} = solve_for_kA;
    mycell = {klabels{:};kcell{:}};
    fixedk = [fixedk,fp{1}.kA,false];
    fixedklabels = {fixedklabels{:},'kA','solve_for_kA'};
    my_utility2(mycell{:},'parcomp',parcomp,'fixedk',fixedk,...
      'fixedklabels',fixedklabels);
    fprintf(['R^2 = %g (normal space ave); \n',...
      'R^2 = %g (log-transformed space ave); \n',...
      'R^2 = %g (normal space all); \n',...
      'R^2 = %g (log-transformed space all);\n'],...
      1-stats.SSEave/stats.SSTave,1-stats.SSElogave/stats.SSTlogave,...
      1-stats.SSEall/stats.SSTall,1-stats.SSElogall/stats.SSTlogall);
    
%     [val,sepsis_exp,sol] = model_scaled_2018_loop(200,100,exd,...
%       k,klabels,'plot_solution',plot_solution,'parcomp',parcomp);
    if eval_stats
      modelfunave = @(kd,Btd) model_scaled_2018_loop(Btd(:,1),Btd(:,2),exd,...
        kd,klabels,'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels);
      modelfunall = @(kd,Btd) model_scaled_2018_loop(Btd(:,1),Btd(:,2),exd,...
        kd,klabels,'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels);
      expparamodelfunave = @(kexpd,Btd) model_scaled_2018_loop(...
        Btd(:,1),Btd(:,2),exd,k.*exp(kexpd-1),...
        klabels,'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels);
      expparamodelfunall = @(kexpd,Btd) model_scaled_2018_loop(...
        Btd(:,1),Btd(:,2),exd,k.*exp(kexpd-1),...
        klabels,'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels);
      logdatamodelfunave = @(kd,Btd) log(...
        model_scaled_2018_loop(Btd(:,1),Btd(:,2),exd,kd,klabels,...
        'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels));
      logdatamodelfunall = @(kd,Btd) log(...
        model_scaled_2018_loop(Btd(:,1),Btd(:,2),exd,kd,klabels,...
        'plot_solution',false,'parcomp',parcomp,...
        'fixedk',fixedk,'fixedklabels',fixedklabels));
      opts = statset('Display','iter','TolFun',1e-10,'MaxIter',0,...
        'UseParallel',true);
      names = {'ave','logave','expave','all','logall','expall'};
      inds = Bstave(:,2) > 0
      clear kc;
      Bst{1} = Bstave; B{1} = Bave; mf{1} = modelfunave; kc{1} = k;
      Bst{2} = Bstave(inds,:); B{2} = log(Bave(inds,:)); kc{2} = k;
      mf{2} = logdatamodelfunave;
      Bst{3} = Bstave; B{3} = Bave; mf{3} = expparamodelfunave;
      kc{3} = ones(size(k));
      Bstall = [exd.Bs.list',exd.t.list']; Ball = exd.B.list';
      indsall = Bstall(:,2) > 0;
      Bst{4} = Bstall; B{4} = Ball; mf{4} = modelfunall; kc{4} = k;
      Bst{5} = Bstall(indsall,:); B{5} = log(Ball(indsall,:));
      mf{5} = logdatamodelfunall; kc{5} = k;
      Bst{6} = Bstall; B{6} = Ball; mf{6} = expparamodelfunall;
      kc{6} = ones(size(k));
%       for mc = 1:numel(mf)
%         mdl{mc} = fitnlm(Bst{mc},B{mc},mf{mc},kc{mc},'Options',opts);
%       end
      modelind = 1;
      if log_transformed, modelind = modelind+1; end
      if use_all, modelind = modelind+3; end
      mdl{modelind} = fitnlm(Bst{modelind},B{modelind},mf{modelind},...
        kc{modelind},'Options',opts);
      coefCI(mdl{modelind})
    end
    
    close(figure(9)); figure(9);
    subplot(2,2,1);
    for sc = 1:numel(sol)
      plot(sol{sc}.solution_time,sol{sc}.bacteria,[mycols(sc),'-']);
      hold on
    end
    plot(exd.t.list,mean([exd.B.ave{:}]).*ones(size(exd.t.list)),'k-',...
      exd.t.list,mean(exd.B.list).*ones(size(exd.t.list)),'k:');
    for sc = 1:numel(sol)
      plot(exd.t.ave{sc},exd.B.ave{sc},[mycols(sc),'*']);
      hold on
    end
%     p1 = plot(sol{1}.solution_time,sol{1}.bacteria,'r-',...
%       sol{2}.solution_time,sol{2}.bacteria,'b-',...
%       sol{3}.solution_time,sol{3}.bacteria,'g-',...
%       sol{4}.solution_time,sol{4}.bacteria,'k-',...
%       exd.t.list,mean([exd.B.ave{:}]).*ones(size(exd.t.list)),'k-',...
%       exd.t.list,mean(exd.B.list).*ones(size(exd.t.list)),'k:',...
%       exd.t.ave{1},exd.B.ave{1},'r*',...
%       exd.t.ave{2},exd.B.ave{2},'b*',...
%       exd.t.ave{3},exd.B.ave{3},'g*',...
%       exd.t.ave{4},exd.B.ave{4},'k*');
    legentry = {};
    for sc = 1:numel(sol)
      legentry{sc} = ['B_0 = ',num2str(exd.Bsource(sc),'%3.2g'),'-Actual'];
    end
    legend(legentry{:},'Average of the average','Average of all');
    xlim([0,100]);
    subplot(2,2,2);
    for sc = 1:numel(sol)
      plot([exd.t.all{sc}{:}],[exd.B.all{sc}{:}],[mycols(sc),'*'],...
        sol{sc}.solution_time,sol{sc}.bacteria,[mycols(sc),'-']);
      hold on
    end
    plot(exd.t.list,mean([exd.B.ave{:}]).*ones(size(exd.t.list)),'k-',...
      exd.t.list,mean(exd.B.list).*ones(size(exd.t.list)),'k:');
%     plot([exd.t.all{1}{:}],[exd.B.all{1}{:}],'r*',...
%       sol{1}.solution_time,sol{1}.bacteria,'r-',...
%       [exd.t.all{2}{:}],[exd.B.all{2}{:}],'b*',...
%       sol{2}.solution_time,sol{2}.bacteria,'b-',...
%       [exd.t.all{3}{:}],[exd.B.all{3}{:}],'g*',...
%       sol{3}.solution_time,sol{3}.bacteria,'g-',...
%       [exd.t.all{4}{:}],[exd.B.all{4}{:}],'k*',...
%       sol{4}.solution_time,sol{4}.bacteria,'k-',...
%       exd.t.list,mean([exd.B.ave{:}]).*ones(size(exd.t.list)),'k-',...
%       exd.t.list,mean(exd.B.list).*ones(size(exd.t.list)),'k:');
%     hold on
    for sc = 1:numel(exd.t.allo)
      plot([exd.t.allo{sc}],[exd.B.allo{sc}],[mycols(sc),'o']);
    end
    xlim([0,100]);
    subplot(2,2,3);
    for sc = 1:numel(exd.t.ave)
      plot(exd.t.ave{sc},log(exd.B.ave{sc}),[mycols(sc),'*'],...
        sol{sc}.solution_time,log(sol{sc}.bacteria),[mycols(sc),'-']);
      hold on
    end
    plot(exd.t.list,mean(log([exd.B.ave{:}])).*ones(size(exd.t.list)),'k-',...
      exd.t.list,mean(log(exd.B.list)).*ones(size(exd.t.list)),'k:');
    ylim([min(log(exd.B.list)),max(log(exd.B.list))]);
    xlim([0,100]);
    subplot(2,2,4);
    for sc = 1:numel(exd.t.all)
      plot([exd.t.all{sc}{:}],log([exd.B.all{sc}{:}]),[mycols(sc),'*'],...
        sol{sc}.solution_time,log(sol{sc}.bacteria),[mycols(sc),'-']);
      hold on
    end
    plot(exd.t.list,mean(log([exd.B.ave{:}])).*ones(size(exd.t.list)),'k-',...
      exd.t.list,mean(log(exd.B.list)).*ones(size(exd.t.list)),'k:');
    for sc = 1:numel(exd.t.allo)
      plot([exd.t.allo{sc}],log([exd.B.allo{sc}]),'ro');
    end
    ylim([min(log(exd.B.list)),max(log(exd.B.list))]);
    xlim([0,100]);
    
    if save_for_debug
      set(gcf,'Units','inches','Position',[1,1,8,6]);
      saveas(gcf,'Bacterial_fits.svg');
      
      close(figure(101)); figure(101);
      set(gcf,'Units','inches','Position',[1,1,12,9]);
      
      fn = {'bacteria','pro_inf','anti_inf','damage','clot','bflux',...
        'aprod','eq3a'};
      if plot_b_terms
        fn = {'bacteria','pro_inf','anti_inf','damage','clot','Blog',...
          'loci','gloi'};
      end

      for Bloadc = 1:numel(sol)
        for fnc = 1:8
          subplot(3,3,fnc);
          plot(sol{Bloadc}.solution_time,sol{Bloadc}.(fn{fnc}),...
            mycols(Bloadc));
          hold on;
          title(fn{fnc});
          xlim([0,200]);
        end
      end
      if plot_b_terms
        for Bloadc = 1:numel(sol)        
          ax = subplot(3,3,Bloadc+5); delete(ax); subplot(3,3,Bloadc+5);
          mystys = {':','-.','--','-'};
          for fnc = 6:8
            plot(sol{Bloadc}.solution_time,sol{Bloadc}.(fn{fnc}),...
              [mycols(Bloadc),mystys{fnc-5}]);
            hold on;
          end
          xlim([0,200]);
        end
        legend('Logistic','Local','Global');
      end
      legentry = {};
      for Bsc = 1:numel(exd.Bsource)
        legentry{Bsc} = ['B_0 = ',num2str(exd.Bsource(Bsc),'%3.2g')];
      end
      subplot(3,3,1); legend(legentry); saveas(gcf,'Dynamics.svg');
      
      close(figure(102)); figure(102);
      set(gcf,'Units','inches','Position',[1,1,8,6]);
      
      ini_ref = 2; nref = 2;
      mycell = {klabels{:};kcell{:}};
      vars = {'k5','nu1','kBl','kA','k1','c1'};
%       vars = {'k1'};
      close(figure(102)); figure(102);
      fn = {'bacteria','pro_inf','anti_inf','damage'};
      for fnc = 1:numel(fn), close(figure(200+fnc)); end
      lower_fact = 1/2;  upper_fact = 10;
      lower_fact = 0;  upper_fact = 2;
      for vc = 1:numel(vars)
        nr = ceil(sqrt(numel(vars))); nc = ceil(numel(vars)/nr);
        inds = ([1:numel(vars)] ~= vc);
        othervarnames = vars(inds);        
        myvar = kcellfull{1}(ismember(klabelcellfull{1},vars{vc}));
        othervar = kcellfull{1}(ismember(klabelcellfull{1},othervarnames));
        inds2 = [];
        for myc = 1:size(mycell,2)
          if ~ismember(vars{vc},mycell{1,myc})
            inds2 = [inds2,myc];
          end
        end
        mycell2 = {mycell{1,inds2};mycell{2,inds2}};
        f = @(var,Bsource_in) my_finder(vars{vc},var,...
          'Bsource_in',Bsource_in,mycell{:},'fixedk',[fixedk,othervar],...
          'fixedklabels',{fixedklabels{:},othervarnames});
        
        figure(103);
        [xg{vc},yg{vc},fg{vc},fgo{vc}] = interface_finder('xl',lower_fact*myvar,...
          'xu',upper_fact*myvar,'yl',0,'yu',5e2,'f',f,'noutputs',2,...
          'ini_ref',ini_ref,'nref',nref);
        figure(102); subplot(nr,nc,vc);
        [ci{vc},hi{vc}] = contourf(xg{vc},yg{vc},fg{vc},[-0.5,0.5,1.5]);
        xlabel(vars{vc}); ylabel('Bsource');
        
        for fnc = 1:numel(fn)
          figure(200+fnc); subplot(nr,nc,vc);
          final.(fn{fnc}) = fg{vc};
          for lc = 1:numel(fgo{1})
            final.(fn{fnc})(lc) = fgo{vc}{lc}.(fn{fnc})(end);
          end
          if range(final.(fn{fnc}),'all') > 0
            contourf(xg{vc},yg{vc},final.(fn{fnc})); colorbar
          end
          %  For debugging
          figure(300+fnc); subplot(nr,nc,vc);
          for lc = 1:numel(fgo{1})
            plot(fgo{vc}{lc}.solution_time,fgo{vc}{lc}.(fn{fnc})); hold on;
          end
        end
      end
      for fnc = 1:numel(fn)
        figure(200+fnc); 
        title(fn{fnc});
      end
      
    end

    if eval_stats
%       mdl{modelind}
      for mc = modelind% = 1:numel(mdl)
        mdl{mc}
        close(figure(9+mc+fig_bump)); figure(9+mc+fig_bump);
        set(gcf,'Units','inches','Position',[1,1,14,6]);
        
        if eval_stats, set(gcf,'Name',names{mc}); end
        
        subplot(2,3,1);
        for sc = 1:numel(sol)
          plot(sol{1}.solution_time,sol{1}.bacteria,[mycols(sc),'-']);
          hold on
        end
        plot(exd.t.list,mean([exd.B.ave{:}]).*ones(size(exd.t.list)),'k-',...
          exd.t.list,mean(exd.B.list).*ones(size(exd.t.list)),'k:');
        legentry = {};
        for ec = 1:numel(exd.t.ave)
          plot(exd.t.ave{ec},exd.B.ave{ec},[mycols(ec),'*']);
          legentry{ec} = ['B_0 = ',num2str(exd.Bsource(ec),'%3.2g')];
        end
        legend(legentry{:},'Average of the average','Average of all');
        
        subplot(2,3,2);
        normplot(mdl{mc}.Residuals.Standardized);
        %  Estimates regarding normality of the residuals, Please cite this as:
        %  �ner, M., & Deveci Kocako�, ?. (2017). JMASM 49: A Compilation of
        %  Some Popular Goodness of Fit Tests for Normal Distribution: Their
        %  Algorithms and MATLAB Codes (MATLAB). Journal of Modern Applied
        %  Statistical Methods, 16(2), 30.
        normalitytest(mdl{mc}.Residuals.Standardized');
        subplot(2,3,3);
        plot(B{mc}-mdl{mc}.Residuals.Raw,B{mc},'.',B{mc},B{mc});
        xlabel('Predicted'); ylabel('Actual');
        subplot(2,3,4);
        plot(B{mc},mdl{mc}.Residuals.Standardized,'+'); grid on
        xlabel('Actual B'); ylabel('Standardized Res');
        subplot(2,3,5);
        plot(Bst{mc}(:,1),mdl{mc}.Residuals.Standardized,'+');
        xlabel('Bsource'); ylabel('Standardized Res');
        subplot(2,3,6);
        plot(Bst{mc}(:,2),mdl{mc}.Residuals.Standardized,'+');
        xlabel('times'); ylabel('Standardized Res');
        
%         specs = {'damage','pro_inf','anti_inf'};
%         for nc = 1:numel(specs)
%           subplot(3,3,6+nc);
%           plot(sol{1}.solution_time,sol{1}.(specs{nc}),'r-',...
%             sol{2}.solution_time,sol{2}.(specs{nc}),'b-',...
%             sol{3}.solution_time,sol{3}.(specs{nc}),'g-',...
%             sol{4}.solution_time,sol{4}.(specs{nc}),'k-');
%           title(specs{nc});
%         end
        if save_for_debug
          saveas(gcf,'Stat_tests.svg');
          standardized_residuals = mdl{mc}.Residuals.Standardized;
          save std_resids standardized_residuals;
        end
      end
      
    end
    
%     figure(9)
      
    if eval_stats & (1 == 0)
      for Bloadc = 1:numel(exd.Bsource)
        figure(10+Bloadc);
        for expc = 1:numel(exd.t.all{Bloadc})
          for tc = 1:numel(exd.t.all{Bloadc}{expc})
            if ~exd.outlier{Bloadc}{expc}(tc)
              ind = find(exd.t.all{Bloadc}{expc}(tc) == exd.t.ave{Bloadc});
              res{Bloadc}{expc}(tc) = abs(exd.B.all{Bloadc}{expc}(tc)-...
                sim{Bloadc}(ind));
            else
              res{Bloadc}{expc}(tc) = nan;
            end
          end
          plot(exd.t.all{Bloadc}{expc},res{Bloadc}{expc}/std(exd.B.list)); hold on;
        end
        legend;
      end

      for Bloadc = 1:numel(exd.Bsource)
        figure(20+Bloadc);
        for expc = 1:numel(exd.t.all{Bloadc})
          for tc = 1:numel(exd.t.all{Bloadc}{expc})
            if ~exd.outlier{Bloadc}{expc}(tc)
              ind = find(exd.t.all{Bloadc}{expc}(tc) == exd.t.ave{Bloadc});
              res{Bloadc}{expc}(tc) = abs(exd.B.all{Bloadc}{expc}(tc)-...
                sim{Bloadc}(ind));
            end
          end
          plot(exd.t.all{Bloadc}{expc},res{Bloadc}{expc}/...
            std([exd.B.all{Bloadc}{:}])); hold on;
        end
        legend;
      end

      for Bloadc = 1:numel(exd.Bsource)
        figure(31);
        for expc = 1:numel(exd.t.all{Bloadc})
          for tc = 1:numel(exd.t.all{Bloadc}{expc})
            if ~exd.outlier{Bloadc}{expc}(tc)
              ind = find(exd.t.all{Bloadc}{expc}(tc) == exd.t.ave{Bloadc});
              res{Bloadc}{expc}(tc) = abs(exd.B.all{Bloadc}{expc}(tc)-...
                mean([exd.B.all{Bloadc}{:}]));
            end
          end
          plot(exd.t.all{Bloadc}{expc},res{Bloadc}{expc}/...
            std([exd.B.all{Bloadc}{:}])); hold on;
        end
        legend;
      end

      for Bloadc = 1:numel(exd.Bsource)
        figure(32);
        for expc = 1:numel(exd.t.all{Bloadc})
          for tc = 1:numel(exd.t.all{Bloadc}{expc})
            if ~exd.outlier{Bloadc}{expc}(tc)
              ind = find(exd.t.all{Bloadc}{expc}(tc) == exd.t.ave{Bloadc});
              res{Bloadc}{expc}(tc) = abs(log(exd.B.all{Bloadc}{expc}(tc))-...
                mean(log([exd.B.all{Bloadc}{:}])));
            end
          end
          plot(exd.t.all{Bloadc}{expc},res{Bloadc}{expc}/...
            std(log([exd.B.all{Bloadc}{:}]))); hold on;
        end
        legend;
      end

      figure(33);
      plot((exd.B.list-mean(exd.B.list))./std(exd.B.list));

      figure(34);
      plot((log(exd.B.list)-log(mean(exd.B.list)))./std(log(exd.B.list)));

    end
    
%       break
%     figure(10)
%     plot(time1,actual_conc1,'r*',time1,sepsis_exp1,'r-',...
%         time2,actual_conc2,'b*',time2,sepsis_exp2,'b-',...
%         time3,actual_conc3,'g*',time3,sepsis_exp3,'g-',...
%         time4,actual_conc4,'k*',time4,sepsis_exp4,'k-');
%     legend(['B_0 = ',num2str(exd.Bsource1,'%3.2g'),'-Actual'],...
%         ['B_0 = ',num2str(exd.Bsource1,'%3.2g'),'-Model'],...
%         ['B_0 = ',num2str(exd.Bsource2,'%3.2g'),'-Actual'],...
%         ['B_0 = ',num2str(exd.Bsource2,'%3.2g'),'-Model'],...
%         ['B_0 = ',num2str(exd.Bsource3,'%3.2g'),'-Actual'],...
%         ['B_0 = ',num2str(exd.Bsource3,'%3.2g'),'-Model'],...
%         ['B_0 = ',num2str(exd.Bsource4,'%3.2g'),'-Actual'],...
%         ['B_0 = ',num2str(exd.Bsource4,'%3.2g'),'-Model']);
%     
%     figure(11)
%     hold on
%     plot(sol1.solution_time, sol1.sepsis,'r','linewidth',2)
%     plot(sol2.solution_time, sol2.sepsis,'b','linewidth',2)
%     plot(sol3.solution_time, sol3.sepsis,'g','linewidth',2)
%     plot(sol4.solution_time, sol4.sepsis,'m','linewidth',2)
%     xlabel('time')
%     ylabel('Bacteria in lumen')
% 
%     figure(12)
%     hold on
%     plot(sol1.solution_time, sol1.permeability,'r','linewidth',2)
%     plot(sol2.solution_time, sol2.permeability,'b','linewidth',2)
%     plot(sol3.solution_time, sol3.permeability,'g','linewidth',2)
%     plot(sol4.solution_time, sol4.permeability,'m','linewidth',2)
%     xlabel('time')
%     ylabel('Permeability')
%     
%     figure(13)
%     hold on
%     plot(sol1.solution_time, sol1.bacteria_blood,'r','linewidth',2)
%     plot(sol2.solution_time, sol2.bacteria_blood,'b','linewidth',2)
%     plot(sol3.solution_time, sol3.bacteria_blood,'g','linewidth',2)
%     plot(sol4.solution_time, sol4.bacteria_blood,'m','linewidth',2)
%     xlabel('time')
%     ylabel('Bacteria in blood')
% 
%     figure(14)
%     hold on
%     plot(sol1.solution_time, sol1.M_cells,'r','linewidth',2)
%     plot(sol2.solution_time, sol2.M_cells,'b','linewidth',2)
%     plot(sol3.solution_time, sol3.M_cells,'g','linewidth',2)
%     plot(sol4.solution_time, sol4.M_cells,'m','linewidth',2)
%     xlabel('time')
%     ylabel('Immune cells')
% 
%     figure(15)
%     hold on
%     plot(sol1.solution_time, sol1.A_cells,'r','linewidth',2)
%     plot(sol2.solution_time, sol2.A_cells,'b','linewidth',2)
%     plot(sol3.solution_time, sol3.A_cells,'g','linewidth',2)
%     plot(sol4.solution_time, sol4.A_cells,'m','linewidth',2)
%     xlabel('time')
%     ylabel('Anti-inflammatory cells')
    
end
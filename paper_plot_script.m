%%  %  Some initial data processing.  Need to run before making plots.
%  Running this script/pressing F5 will process all data up until the first
%  "return" below.  Data from previous figures may be needed for subsequent
%  figures, so run each figure from start to finish to see them all.
%%  This folder ignores 1.94e? data
%  The data has uncertainties associated with it making it unsuitable for
%  use in the calibrating data set
outlier_matrix_194 = [4,1,1;4,1,2;4,1,3;4,1,4;4,1,5;4,1,6;...
  4,2,1;4,2,2;4,2,3;4,2,4; 4,3,1;4,3,2;4,4,1];
plot_all_data(outlier_matrix_194);

%%  Optimization and integration toolboxes
%  Location of paths (as of 2019/05/09)
addpath ./optimization/
addpath ./Finalized_Code/
addpath ./Matlab_tools/

%%  This is where we save figures
mydir = '.\';

%%  Whether or not to save the pics
save_pics = false;

%%  My favorite parameter set as of 4/26/2021
% Healthy < Bsource = 128.107 < Aseptic < Bsource = 505 < Septic.
% R^2 = -1.24142 (normal space ave); 
% R^2 = 0.567363 (log-transformed space ave); 
% R^2 = -0.224484 (normal space all); 
% R^2 = 0.70173 (log-transformed space all);
%  Normal 10/10
param_names = {'B_infy','f','kD1','Bmf','clot_capacity','nu4','k1'};
param_vals = [0.736964699290478
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
fdebug_script_paper;

%%  Find/deal with kA and fixed (run only once per above parameter values)
%  We pretty much always solve for kA so that way anti-inflammatory
%  inhibition is 75%:  1/(1+kA*Amax) = 0.75.  Thus we rarely know what kA
%  is as it is actually changing to make sure 1/(1+kA*Amax) is always 0.75
%  = inhi_perc (which theoretically we can change if we want, but we
%  haven't explored other values of inhi_perc at this time).
[my_quantity,stats,sol] = fdebugvarg(param_vals,'plot_sol',false);
param_names = {param_names{:},'kA'};
param_vals(end+1) = sol{3}.struc.fp.kA;
for ac = 1:numel(addlargs)
  if isequal(addlargs{ac},'solve_for_kA')
    addlargs{ac+1} = false;
  end
end

%  Get names of all parameters (and their descriptions)
for aac = 1:numel(addlargs)
  if isequal(addlargs{aac},'fixedk')
    param_vals = [param_vals,addlargs{aac+1}];
  elseif isequal(addlargs{aac},'fixedklabels')
    param_names = {param_names{:},addlargs{aac+1}{:}}
  end
end

fdebug_script_paper;

%%  Read in the plot-digitizer data
%  In case I want to use this on Linux at some point
if ispc, my_slash = '\'; else, my_slash = '/'; end

%  Reads the data from the csv file (converted from the students' original
%  xml file using plot digitizer, 1/22/2020) starting from row 7 (before
%  that is just plot digitizer header junk)
% M = csvread(['..',my_slash,'..',my_slash,'Rat Mortality',my_slash,...
%   'Data_points.csv'],6);
M = csvread('Rat_Mortality_Data_points.csv',6);

%  Get "bacteria implanted" level estimates.   Also, renormalize.  Still
%  not sure how best to renormalize these (there was discussion of CFUs,
%  bacteria numbers, grams, ccs).  I divide by 1e6 for now.
bacteria_levels_init = M(:,2)./1e6;
%  Corresponding "mortality" data (read it off the chart)
mortality_data_init = [24 24 24 24 48 24 24 24 48 48 48 48 72 48 ...
  72 72 72 72 24 96 96 48 48 48 48 96 96];
outcome = mortality_data_init > 24;
possible_hrs = [24,48,72,96];

%  The indices in the picture are weird.  I list them here for later ease
%  with trying out different runs.
indices = [1:19 21:26 29 30];

%  For easiest access, I went ahead and created 30 component vectors
%  placing 0s in a component when there is no corresponding rat.  That is
%  bacteria_levels(20) = bacteria_levels(27) = bacteria_levels(28) = 0
%  since there is no R-20, R-27, or R-28 in the plot.
bacteria_levels = zeros(30,1);
bacteria_levels(indices) = bacteria_levels_init;
mortality_data(indices) = mortality_data_init;

%%  Statistical analyses for mortality data
%  anova (and ttest2) compares means assuming same variances (and normal
%  distributions)
%  We first test the equal variances assumption:
%    The p-value tells us about the probability that one of the variances
%    is, in fact, significantly different than the others.  With 0.68, it
%    seems like we cannot assert that we have different variances.
[p,stats] = vartestn(bacteria_levels_init,mortality_data_init);

%  We now test the normality assumption following Devore's suggestion of
%  grouping them all into one group after shifting them by their mean
%    The normality tests are 10/10 suggesting normality may be a reasonable
%    assumption.  The normality plot looks fine but I don't have much
%    experience in looking at such normality plots by eye (hence the
%    reliance on the normality tests)
close(figure(2+10)); figure(2+10);
possible_hrs = 24:24:96;
mylist = [];
for phc = 1:numel(possible_hrs)
  inds = mortality_data_init == possible_hrs(phc);
  groups(phc).raw = bacteria_levels_init(inds);
  groups(phc).shif = groups(phc).raw-mean(groups(phc).raw);
  mylist = [mylist;groups(phc).shif];
end
normalitytest(mylist')
normplot(mylist)

%  We now proceed with anova using Matlab's anova1
%    The extemely low p-value suggests at least one of them has a different
%    mean compared to the others.  Based on our past ttest2 applications
%    (which compares just two samples, we already expect that the OMT 24h
%    is the primary offender.
[p,anovatab,stats] = anova1(bacteria_levels_init,mortality_data_init);
c = multcompare(stats, 'alpha', 0.05)

%  Make a table
for ic = 1:numel(bacteria_levels_init)
  fprintf('%2d %3d %2d\n',ic,round(bacteria_levels_init(ic)),...
    mortality_data_init(ic));
end

%%  We now try other possible combinations
%    The first one is the >24 hr group.  Anova yields a p-value of 0.2509
%    suggesting we can't really say that their means are significantly
%    different.
inds = mortality_data_init ~= 24;
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))
%  Honestly, I don't see much point in trying the others, but we go ahead
%  anyway.  We expect all combinations involving OMT 24h will result in a
%  p-value < 0.05 suggesting that one of the groups (OMT 24h) has a
%  significantly different mean than the others.  This seems to be true.
inds = mortality_data_init ~= 48;
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))
inds = mortality_data_init ~= 72;
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))
inds = mortality_data_init ~= 96;
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))

inds = (mortality_data_init ~= 48) & (mortality_data_init ~= 72);
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))
inds = (mortality_data_init ~= 48) & (mortality_data_init ~= 96);
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))
inds = (mortality_data_init ~= 96) & (mortality_data_init ~= 72);
p = anova1(bacteria_levels_init(inds),mortality_data_init(inds))

return

%%  Figure 1
%  In version 2, we first plot the data
close(figure(1)); figure(1);
[exd] = get_data(outlier_matrix_194,'td',0);
mycols = 'gbrk';
xbuff = 0.11; xspan = 0.37; xgap = 0.5-(xbuff+xspan);
ystart = 0.2; yspan = 0.75;
ax1 = axes('Position',[xbuff,ystart,xspan,yspan]);
% subplot(1,2,2);
%  Matlab lies about how it renders figures.  If you say I want a 5x6
%  figure, it somehow converts that into pixels, inaccurately, and uses
%  that.  This factor theoretically fixes that issue and is specific to
%  each computer.
compfact = 1;%1.167;%6/5.25;
linw = compfact*1;  %  LineWidth (in pts)
fons = compfact*12;  %  Fontsize (in pts)
figh = compfact*3;
figw = 2*compfact*3;
fonnam = 'Times New Roman';
fonnam = 'Arial';

for Bloadc = 1:numel(exd.Bsource)
  t = exd.t.ave{Bloadc};
  B = exd.B.lave{Bloadc};
  Bmin = exd.B.min{Bloadc};
  Bmax = exd.B.max{Bloadc};
%   han = errorbar(t,B,B-Bmin,Bmax-B);
%   set(han,'Color',mycols(Bloadc));
  plot(t,B,mycols(Bloadc));
  hold on
end

for Bloadc = 1:numel(exd.Bsource)
  t = [exd.t.all{Bloadc}{:}];
  B = [exd.B.all{Bloadc}{:}];
  plot(t,B,[mycols(Bloadc),'.'],'MarkerSize',15);
end
set(gca,'YScale','log');

set(gca,'YScale','log','FontSize',fons,'FontName',fonnam,'LineWidth',linw);
% set(gca,'YScale','linear','FontSize',fons,'FontName',fonnam,'LineWidth',linw);
xlabel('Time (hours)');
ylabel('Bacteria levels (10^6 bacteria)');
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
set(gcf,'Units','inches','Position',[1 1 figw figh]);
name = 'Figure1';
a = xlim; b = ylim; text(a(1)+0.064*diff(a),b(2)-0.3*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);

ax2 = axes('Position',[2*xbuff+xspan+xgap,ystart,xspan,yspan]);

%  Makes Box and Whiskers plot
boxchart(mortality_data_init./24,bacteria_levels_init);
hold on;
plot(mortality_data_init./24,bacteria_levels_init,'k.','LineWidth',linw,...
  'MarkerSize',15);

set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw,...
  'Position',[2*xbuff+xspan+xgap,ystart,xspan,yspan]);
xlabel('Observed mortality time (hours)');
ylabel('Injected bacteria (10^6 bacteria)');
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
set(gca,'XTick',[1,2,3,4],'XTickLabels',{'24','48','72','96'});
xlim([0.5,4.5]);
title(''); box on;
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.035*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);

% %  Instead makes a "Tukey's procedure" plot including 95% confidence
% %  intervals based on the entire data set.  Matlab is not entirely clear
% %  with how it produces the "Multiple comparison of means" plots so I have
% %  elected to not use this (4/12/21).  To help understand more, for 4
% %  groups, there are 6 different possible comparisons and confidence
% %  intervals that can be calculated based on those comparisons.  One way to
% %  plot this is for each group, average its 3 confidence
% %  intervals/corresponding standard deviations to make a confidence
% %  interval plot.  There may be other ways to combine the standard
% %  deviations including taking into account 
% [p,anovatab,stats] = anova1(bacteria_levels_init,mortality_data_init,'off');
% figure(2);
% [c,m] = multcompare(stats,'Display','off');
% errorbar(m(:,1),1:4,m(:,2),'horizontal','o');
% hold on;
% plot(bacteria_levels_init,mortality_data_init./24,'kx','LineWidth',linw);
% 
% set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw,...
%   'Position',[2*xbuff+xspan+xgap,ystart,xspan,yspan]);
% ylabel('Observed mortality time (hours)');
% xlabel('Injected bacteria (10^6 bacteria)');
% tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
% tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
% set(gca,'YTick',[1,2,3,4],'YTickLabels',{'24','48','72','96'});
% ylim([0.5,4.5]);
% title(''); box on;
% a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.035*diff(b),'B',...
%   'FontSize',fons,'FontName',fonnam);

if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 2
%  This is the process diagram for the model.  Ally and Amy originally made
%  this in Microsoft publisher but I switched it over to powerpoint
%  (because that's what I am used to, no other reason).  It can be found in
%  Figure/Some_pics.pptx.  One can use powerpoint's "save as" capability to
%  save this as, as is currently the case, Diagram.emf.  emf seems to be
%  the best format for Word.

%%  Figure 3
try
  load('std_resids');
catch
  fdebug(param_vals);
  load('std_resids');
end
%  Plot of bacterial populations vs time plus experimental measurements
close(figure(3)); figure(3);
xbuff = 0.11; xspan = 0.37; xgap = 0.5-(xbuff+xspan);
ystart = 0.2; yspan = 0.75;
ax1 = axes('Position',[xbuff,ystart,xspan,yspan]);
% subplot(1,2,2);
%  Matlab lies about how it renders figures.  If you say I want a 5x6
%  figure, it somehow converts that into pixels, inaccurately, and uses
%  that.  This factor theoretically fixes that issue and is specific to
%  each computer.
compfact = 1;%1.167;%6/5.25;
linw = compfact*1;  %  LineWidth (in pts)
fons = compfact*12;  %  Fontsize (in pts)
figh = compfact*3;
figw = 2*compfact*3;
fonnam = 'Times New Roman';
fonnam = 'Arial';

% 'plot_what' can be 'all' or 'ave'
plot_all_points_2018('k',param_vals,'klabels',param_names,addlargs{:},...
  'plot_what','all');
set(gca,'YScale','log','FontSize',fons,'FontName',fonnam,'LineWidth',linw);
% set(gca,'YScale','linear','FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
set(gcf,'Units','inches','Position',[1 1 figw figh]);
name = 'all_bacteria_dynamics';
name = 'Figure3';
legend off;
ylim([0.1,2*1000])
ylim([0.1,800])
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.3*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);

ax2 = axes('Position',[2*xbuff+xspan+xgap,ystart,xspan,yspan]);

hands = normplot(standardized_residuals);
set(hands(1),'Marker','.','MarkerSize',15);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw,...
  'Position',[2*xbuff+xspan+xgap,ystart,xspan,yspan]);
xlabel('Standardized Residuals');  ylabel('Normal Probability');
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
% xlim([-2.5,2.5]);
% set(gca,'XTick',[-2:2]);
title(''); box on;
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.035*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);

if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 4
%  Matlab lies about how it renders figures.  If you say I want a 5x6
%  figure, it somehow converts that into pixels, inaccurately, and uses
%  that.  This factor theoretically fixes that issue and is specific to
%  each computer.
ybump = 0.05;
rxbump = 0;

k1s = [1.4,1.3,1.2];

compfact = 1;%6/5.25;

figh = compfact*5;  %  Figure height (in inches)
figw = compfact*6;  %  Figure width (in inches)

bini = 128;

%  At one point this was supposed to produce all three cases.  Right now I
%  believe we have three healthy cases...kcurr should be adjusted for that.
fignum = 4; close(figure(fignum)); figure(fignum);
myc = 'rbg';
for pc = 1:numel(k1s)
  inds = find(~ismember(param_vals,'k1'));
  [sepsis_exp,sol] = model_scaled_2018('k',param_vals(inds),...
    'klabels',param_names(inds),'Bsource_in',bini,...
    'k1',k1s(pc),'fn',fignum,'lc',myc(pc),'ls','-',...
    'cf',compfact,addlargs{:});
end
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
set(gcf,'Units','inches','Position',[1 1 compfact*figw compfact*figh]);
name = 'typical_dynamics';
name = 'Figure4';

%  Adjust the axes for more space
myaxs = get(gcf,'Children');

%  Top left
axes(myaxs(4));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
ylim([0,600]);
set(gca,'YScale','log'); ylim([0.01,600]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);

%  Top right
axes(myaxs(3));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
a = xlim; b = ylim;
xlim(a); ylim(b);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom left
axes(myaxs(2));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ylim([0.225,0.475]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'C',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom right
axes(myaxs(1));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ylim([0,600]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'D',...
  'FontSize',fons,'FontName',fonnam);

set(gcf,'Renderer','painters');
if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 5
%  Estimate damage levels at the 24,48,72,96 etc--all
addpath('..\Finalized_Code\');
bscale = 1;%0.85;
clear sols;
for blic = 1:numel(bacteria_levels_init)
  [cats(blic),tcs(blic),maxdams(blic),sols(blic),dampts(blic)] = ...
    mortality_estimate(...
    'k',param_vals,'klabels',param_names,...
    'Bsource_in',bscale*bacteria_levels_init(blic),...
    'plot_solution',false,'maxt',200);
%   fc1 = @(c1,k1) my_finder('c1',c1,'k1',k1,'Bsource_in',Bsource_in,mycell{:});
%   sah(blic) = my_finder('Bsource_in',bscale*bacteria_levels_init(blic));
end

for blic = 1:numel(bacteria_levels_init)
  sols(blic).cumeps = cumtrapz(sols(blic).solution_time,sols(blic).damage);
  sols(blic).cumb = cumtrapz(sols(blic).solution_time,sols(blic).bacteria);
  sols(blic).cumm = cumtrapz(sols(blic).solution_time,sols(blic).pro_inf);
  sols(blic).cuma = cumtrapz(sols(blic).solution_time,sols(blic).anti_inf);
end

damptsvec.d24 = [dampts.d24];
damptsvec.d48 = [dampts.d48];
damptsvec.d72 = [dampts.d72];
damptsvec.d96 = [dampts.d96];

close(figure(2)); figure(2);
plot(bacteria_levels_init,damptsvec.d24,...
  bacteria_levels_init,damptsvec.d48,...
  bacteria_levels_init,damptsvec.d72,...
  bacteria_levels_init,damptsvec.d96)

close(figure(3)); figure(3);
plot(damptsvec.d24); hold on;
plot(damptsvec.d48); plot(damptsvec.d72); plot(damptsvec.d96);
fprintf('Done\n');

%%  Figure 5 actual plotting
%  Try the easier statistic assessments that we think make sense
%  Confusion matrix
fignum = 5;
ecrit_vals = linspace(0,500,10001);
bscrit_vals = linspace(0,600,2001);

inds = 1:numel(outcome);
outliers = [];
% outliers = [7,18];
inds = setdiff(inds,outliers);
%  Assume #18 is an outlier
expd = outcome(inds);
simdam = damptsvec.d24(inds);
bli = bacteria_levels_init(inds);
mdi = mortality_data_init(inds);
% sahtmp = sah(inds);
solstmp = sols(inds);

%  Two options, use the model/epsilon_critical or just use the initial
%  bacterial load for prediction
mypred = 'ecrit';
switch mypred
  case 'bscrit'
    crit_vals = bscrit_vals;
  case 'ecrit'
    crit_vals = ecrit_vals;
end

for cvc = 1:numel(crit_vals)
  
  switch mypred
    case 'ecrit'
      sim = simdam<crit_vals(cvc);
    case 'bscrit'
      sim = bli'<crit_vals(cvc); 
  end
  
  confusion_matrix{cvc} = [sum(expd == true & sim == true),...
    sum(expd == true & sim == false);
    sum(expd == false & sim == true),...
    sum(expd == false & sim == false)]';
  confusion_matrix{cvc} = confusion_matrix{cvc}./sum(confusion_matrix{cvc}(:));
  accuracy(cvc) = sum(diag(confusion_matrix{cvc}));
  
  sse = sum((sim-expd).^2);
  sst = sum((sim-mean(sim)).^2);
  R2(cvc) = 1-sse/sst;
end

tmp = max(accuracy);
inds = accuracy == tmp;
max_crit_val = ceil(mean(crit_vals(inds)));
sim = simdam < max_crit_val;

explivefrac = sum(expd)/numel(expd);
expdeadfrac = 1-explivefrac;
simlivefrac = sum(sim)/numel(sim);
simdeadfrac = 1-simlivefrac;

chi2 = (explivefrac-simlivefrac)^2/simlivefrac+...
  (expdeadfrac-simdeadfrac)^2/simdeadfrac;

%  AUC-ROC type stuff.  Throughout, positive means the animals lived past
%  24 h.
%  Easier access
all = reshape([confusion_matrix{:}],2,2,[]);
%  True positives (predicted positives are actually positives)
TP = squeeze(all(1,1,:));
%  False negatives (predicted negatives are actually positives)
FN = squeeze(all(2,1,:));
%  False positives (predicted positives are actually negatives)
FP = squeeze(all(1,2,:));
%  True negatives (predicted negatives are actually negatives)
TN = squeeze(all(2,2,:));

%  TPR/recall/sensitivity (fraction of actual trues that are correctly
%  predicted)
TPR = TP./(TP+FN);
%  FPR = 1-sensitivity (sensitivity = fraction of actual falses that are
%  correctly predicted, FPR = fraction of actual falses that are
%  incorrectly predicted)
FPR = FP./(TN+FP);

[chim,chiv] = chi2stat(numel(sim)-1);
%  "Degrees of freedom" = 2 and corresponds to the number of levels to be
%  tested...2 here.  Note this assumes # of unknown/determined parameters
%  is 0 (technically true, I think, because we didn't do fitting using this
%  data, we used the other data), that "n_i*p_i" = frequency of the ith
%  event occuring > 5 (i.e. "large n"), and that initial bacterial levels
%  don't matter (which they do).  The conclusion from this test, however,
%  is that our hypothesized "probability distribution" cannot be rejected.
%  The "distribution" assumed lack of dependence on initial bacterial
%  levels and is produced using our model plus the critical damage level
%  from above.
p = chi2cdf(chi2,1);
fprintf(['Chi2 test (no initial bacterial dependence assumed) suggests\n',...
  'p-value for simulation hypothesized distribution rejection is %g\n'],...
  1-p);

%  Cumulative distribution functions (just sum up all successes
[blisorted,inds] = sort(bli);
simcdf = cumsum(sim(inds)); simcdf = simcdf/simcdf(end);
expcdf = cumsum(expd(inds)); expcdf = expcdf/expcdf(end);
close(figure(1)); figure(1);
plot(blisorted,simcdf,blisorted,expcdf);
ai = abs(simcdf-expcdf);
max(ai)
[h,p,ks2stat] = kstest2(bacteria_levels_init(sim==1),...
  bacteria_levels_init(expd==1));
xlabel('Initial bacterial levels');
ylabel('Corresponding cumulative distribution function');
legend('Simulation','Experiment');
% bi = abs(simcdf(2:end)-expcdf(1:end-1))

close(figure(fignum)); figure(fignum);
name = ['Figure',num2str(fignum)];
figh = compfact*3;  %  Figure height (in inches)
figw = compfact*6;  %  Figure width (in inches)
set(gcf,'Units','inches','Position',[1,1,figw,figh]);
vertbump = 0.175;
subplot(1,2,2);
ax1 = get(gca);
ax1h = gca;

subplot(1,2,1);
ax2 = get(gca);
ax2h = gca;

plot(crit_vals,accuracy,'b');
switch mypred
  case 'ecrit'
    xlim([40,150])
%     xlim([0,500])
  case 'bscrit'
    xlim([0,600]);
end
% plot(crit_vals,accuracy,'b',crit_vals,R2,'r');
% ylim([0.5,1])
xlabel('Critical damage level (\epsilon_{crit})');
ylabel('Model prediction accuracy');
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.05*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ys = get(gca,'YTick');
for ysc = 1:numel(ys)
  hold on
  plot(xlim,ys(ysc)*[1 1],'k:');
end

subplot(1,2,2);
plot(FPR,TPR,'.-r','LineWidth',linw);
xlabel('False Positive Rate'); ylabel('True Positive Rate');
% title(sprintf('Area Under Curve = %3.4g',trapz(FPR,TPR)));
fprintf('Area Under Curve = %3.4g',trapz(FPR,TPR));
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.05*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ys = get(gca,'YTick');

% tmppos = ax1.Position;
% set(ax1h,'Position',[tmppos(1),tmppos(2),tmppos(3),tmppos(4)+vertbump]);
% set(ax1h,'YTick',0:0.2:1);
% tmppos = ax2.Position;
% set(ax2h,'Position',[tmppos(1),tmppos(2)+0.75*vertbump,tmppos(3),tmppos(4)-0.5*vertbump]);
axis equal
set(gcf,'renderer','painters');

if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 6
fignum = 6;
close(figure(fignum)); figure(fignum);
set(gcf,'Units','inches','Position',[1,1,6,5]);
xs = 1:numel(bli);
mycs = 'grb';
myst = {'-','-'};{':','--'};
%  Glean kA from the first simulation (assumes kA same for all, which it
%  should be once we have chosen our "optimal" parameter values
tmp = sols(1).struc.extdata.varargin{3};
kA = tmp.kA;
maxt = 96;
simtmp = (sim == true)+1;
for xsc = 1:numel(xs)
  myc = mycs(simtmp(xsc)+1);
  subplot(2,2,1);
  solstmp(xsc).extdata = sols(xsc).struc.extdata;
  solstmp(xsc).y = [solstmp(xsc).bacteria,solstmp(xsc).pro_inf,...
    solstmp(xsc).anti_inf,solstmp(xsc).damage]';
  solstmp(xsc).x = solstmp(xsc).solution_time;
  class = my_classifier(solstmp(xsc))
  plot(solstmp(xsc).solution_time,solstmp(xsc).bacteria,[myc,myst{class}]);
  xlim([0,maxt]); set(gca,'XTick',linspace(0,maxt,5)); grid on; hold on
  ylabel('Bacterial Levels (10^6 bacteria)'); xlabel('Time (hrs)');

  subplot(2,2,2);
  plot(solstmp(xsc).solution_time,solstmp(xsc).pro_inf,myc);
  xlim([0,maxt]); set(gca,'XTick',linspace(0,maxt,5)); grid on; hold on
  ylabel('Pro-inflammatory Levels'); xlabel('Time (hrs)');
  
  subplot(2,2,3);
  plot(solstmp(xsc).solution_time,solstmp(xsc).anti_inf,myc);
  xlim([0,maxt]); set(gca,'XTick',linspace(0,maxt,5)); grid on; hold on
  ylabel('Anti-inflammatory Levels'); xlabel('Time (hrs)');
  
  subplot(2,2,4);
  plot(solstmp(xsc).solution_time,solstmp(xsc).damage,myc);
  xlim([0,maxt]); set(gca,'XTick',linspace(0,maxt,5)); grid on; hold on
  ylabel('Damage Levels'); xlabel('Time (hrs)');
  
end

%  Adjust the axes for more space
myaxs = get(gcf,'Children');
ybump = 0.04;

%  Top left
axes(myaxs(4));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
ylim([0,550]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);

%  Top right
axes(myaxs(3));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
a = xlim; b = ylim;
xlim(a); ylim(b);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom left
axes(myaxs(2));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ylim([0.225,0.425]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'C',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom right
axes(myaxs(1));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
ylim([0,600]);
hold on;
plot(xlim,max_crit_val*[1 1],'m--','LineWidth',2*linw);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'D',...
  'FontSize',fons,'FontName',fonnam);

name = ['Figure',num2str(fignum)];
if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 7
%  Compare rats 11 and 21 "simulations" on same plot
%  You may wish to adjust parameter values to see if you can get the
%  "mortality results" that you might want:
fignum = 7;
%  Rat 11 is OMT 48, Rat 21 is OMT 96
rinds = [11,21];  %rat number
%  With rescaling of initial bacterial levels
crit_val = 20;
% %  Without rescaling of initial bacterial levels
bcctog = 1;
crit_val = 36;
%  Fit k1:  1.27176887084432
k1s = [1.27176887084432,1.2];
nu1s = [0.08,0.08];
lcs = 'rbgk';

close(figure(fignum)); figure(fignum);
clear cats tcs maxdams sols dampts;

name = 'reproduction_of_experimental_variability';
name = 'Figure7';

mycs = 'rg';

for ric = 1:numel(rinds)
%   [sepsis_exp,sol,rp,op,fp,strp] = model_scaled_2018(...
%     'k',param_vals,'klabels',param_names,...
%     'Bsource_in',bcctog*bacteria_levels(rinds(ric)),'c1',c1s(ric),'k1',k1s(ric),...
%     'fn',fignum,'lc',lcs(ric),'ls','-','tspan',[0,96]);
  inds = find(~ismember(param_vals,'k1'));
  [cats(ric),tcs(ric),maxdams(ric),sols(ric),dampts(ric)] = ...
    mortality_estimate(...
    'k',param_vals(inds),'klabels',param_names(inds),...
    'Bsource_in',bacteria_levels(rinds(ric)),'nu1',nu1s(ric),...
    'k1',k1s(ric),'plot_solution',true,'crit_val',crit_val,...
    'fignum',fignum,'lc',mycs(ric),'solve_for_kA',true);
  bcctog*bacteria_levels(rinds(ric))
end

%  Adjust the axes for more space
myaxs = get(gcf,'Children');
ybump = 0.05;
rxbump = 0;

% legend(num2str(rinds'));
%  Top left
axes(myaxs(4));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
% ylim([0,375]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'A',...
  'FontSize',fons,'FontName',fonnam);

%  Top right
axes(myaxs(3));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2)-ybump,tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
xlabel('');
set(gca,'XTickLabel',{});
a = xlim; b = ylim;
xlim(a); ylim(b);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'B',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom left
axes(myaxs(2));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
% ylim([0.225,0.475]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'C',...
  'FontSize',fons,'FontName',fonnam);

%  Bottom right
axes(myaxs(1));
tmppos = get(gca,'Position');
set(gca,'Position',[tmppos(1),tmppos(2),tmppos(3)+rxbump,...
  tmppos(4)+ybump]);
set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
tmp = get(gca,'YLabel'); tmp.FontSize = fons; set(gca,'YLabel',tmp);
tmp = get(gca,'XLabel'); tmp.FontSize = fons; set(gca,'XLabel',tmp);
% ylim([0,175]);
a = xlim; b = ylim; text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),'D',...
  'FontSize',fons,'FontName',fonnam);

set(gcf,'Units','inches','Position',[1 1 6 5]);

for j = 1:4
  axes(myaxs(j));
  for ric = 1:numel(rinds)
%     plot(mortality_data(rinds(ric))*[1 1],ylim,[lcs(ric),':']);
    plot([24;48;72;96]*[1 1],ylim,'k:');
    set(gca,'XTick',[0,24,48,72,96]); xlim([0,96]);
  end
  if j == 1
  plot([0,96],bcctog*max_crit_val*[1 1],'k:');
  tmp = get(gca,'YTick');
  tmp2 = get(gca,'YTickLabel');
  set(gca,'YTick',sort([tmp,max_crit_val]),'YTickLabel',...
    {tmp2{1},'\epsilon_{crit}',tmp2{2:end}});
%   plot([0,96],bcctog*25*[1 1],'k:');
  end
end

set(gcf,'Renderer','painters');
if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end

%%  Figure 8:  calculations for the grids
%  See below sections for prettier renderings of the contours.  For
%  prettier contours increase either the initial level of refinement
%  (ini_ref) or the total number of refinements made (nref).  For the paper
%  ini_ref = 2, nref = 8;
close(figure(1001)); figure(1001);
iniref = 2;
nref = 8;
ns.k1 = [1,1.5];
ns.kBl = [0.4,0.7];
ns.nu1 = [0,0.16];
nBs.kBl = [0,550];
nBs.nu1 = nBs.kBl;
nBs.k1 = nBs.kBl;

cont_names = {'k1','kBl','nu1'};

for pc = 1:numel(cont_names)
  pn = cont_names{pc};
  inds = find(~ismember(param_names,pn));
  kcell = num2cell(param_vals(inds));
  mycell = {param_names{inds};kcell{:}};
  f.(cont_names{pc}) = @(p,Bsource_in) my_finder(pn,p,...
    'Bsource_in',Bsource_in,mycell{:});
  [pvsBx.(pn),pvsBy.(pn),pvsBz.(pn),pvsBo.(pn)] = interface_finder(...
    'xl',ns.(pn)(1),'xu',ns.(pn)(2),'yl',nBs.(pn)(1),'yu',nBs.(pn)(2),...
    'f',f.(pn),'ini_ref',iniref,'nref',nref,'noutputs',2);
end

%%  Figure 8:  Makes plots
%  -alternate make a plot
zoomin = false;
fignum = 8; close(figure(fignum)); figure(fignum);
tmppos = get(gcf,'Position');
set(gcf,'Units','inches','Position',[1,1,3,9]);
% set(gcf,'Units','Inches','Position',1+[0 0 compfact*figw/2 ...
%   compfact*3*figh/2]);
mycols = 'cbr';

for cnc = 1:numel(cont_names)
  pn = cont_names{cnc};
  xs = pvsBx.(pn)(:,1); ys = pvsBy.(pn)(1,:);
  xspad = [xs(1)-(xs(2)-xs(1));xs;xs(end)+(xs(2)-xs(1))];
  yspad = [ys(1)-(ys(2)-ys(1)),ys,ys(end)+(ys(2)-ys(1))];
  cont{cnc}{1} = contourc(xspad,yspad,padarray(pvsBz.(pn),[1 1],1)',[-0.5,0.5]);
  cont{cnc}{2} = contourc(xspad,yspad,padarray(pvsBz.(pn)-2*(pvsBz.(pn) == 2),...
    [1 1],-0.5)',[0.5,1]);
  cont{cnc}{3} = contourc(xspad,yspad,padarray(pvsBz.(pn).*(pvsBz.(pn) == 2),...
    [1 1],0)',[1.5,2.5]);
end

if ~zoomin, rows = 1; else rows = 2; end
for cc1 = 1:numel(cont_names)
  ax{cc1} = subplot(3,1,cc1);
  for cc2 = 1:numel(cont{cc1})
    [~,cx,cy] = extract_contours(cont{cc1}{cc2});
    tmpcell = {};
    for cc = 1:numel(cx)
      tmpcell{3*cc-2} = cx{cc};
      tmpcell{3*cc-1} = cy{cc};
      tmpcell{3*cc} = mycols(cc2);
    end
    fill(tmpcell{:});
    hold on;
  end
end

%%  Figure 8:  Decorate it
%  ybump same as for figure 4
%  xbump changed
%  Bsources hard-wired in for now
Bsources = [128,248,505];
rxbump = -0.04;
lxbump = 0.08;
ybump = 0*0.04;
mylets = 'ABC';

for pc = 1:numel(cont_names)
  pn = cont_names{pc};
  axes(ax{pc});  % Top left
  % xlabel('Bacterial Growth Rate ({\it{k_1}})');
%   if pc == 1
    ylabel('{\it{B_{source}}} (\times 10^6 bacteria/cc)');
%   end
  % set(gca,'XTickLabel',{});
  set(gca,'FontSize',fons,'FontName',fonnam,'LineWidth',linw);
  
  tmppos = get(gca,'Position');
  set(gca,'Position',[tmppos(1)+lxbump,tmppos(2)-ybump,tmppos(3)+rxbump,...
    tmppos(4)+ybump]);
  hold on
  a = ns.(pn); b = nBs.(pn);
  plot([0,2],Bsources'*[1 1],'k:','LineWidth',linw);
  text(a(1)+0.014*diff(a),b(2)-0.06*diff(b),mylets(pc),...
    'FontSize',fons,'FontName',fonnam);
  xlim(a); ylim(b);
  if ~zoomin
    if pc == 2
      xlabel('Local immune response strength ({\it{k_2}})');
%       set(gca,'YTickLabels','');
    elseif pc == 3
      xlabel('Pro-inflammatory activation rate ({\it{\nu_1}})');
%       set(gca,'YTickLabels','');
    elseif pc == 1
      xlabel('Bacterial growth rate ({\it{k_1}})');
    end
  end
end

%  Put x's corresponding to figure 4 runs.
%  ** Should match Figure 2 **
k1s = [1.4,1.3,1.2];
axes(ax{1})
hold on;
plot(k1s,128,'kd','MarkerSize',10,'LineWidth',1.5*linw,...
  'MarkerFaceColor','y');

%  Put o's corresponding to figure 
k1s = [1.27,1.2];
plot(k1s,bacteria_levels(rinds(ric)),'ko','MarkerSize',10,...
  'LineWidth',1.5*linw,'MarkerFaceColor','y');

name = 'Figure8';
set(gcf,'Renderer','painters');
if save_pics
  saveas(gcf,[mydir,name,'.fig']);
  saveas(gcf,[mydir,name,'.emf']);
  saveas(gcf,[mydir,name,'.eps']);
end
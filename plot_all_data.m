function plot_all_data(varargin)
%  Major differences vs 2018:
%  In 2019:
%    1.28e8 B-Experiment 3 omitted
%    5.505e8 B-Experiment 1 4h-190 observation added
%              Experiment 2 4h-800 observation added
%    1.94e9 B-Experiment 1 8h-800 observation added
if round(numel(varargin)/2) ~= numel(varargin)/2
  outlier_matrix = varargin{1}; 
  varargin = {varargin{2:end}};
else
  outlier_matrix = []; 
end

[exd] = get_data(outlier_matrix,varargin{:});

close(figure(2)); figure(2); hold on
close(figure(1)); figure(1); hold on
mycols = 'rgbc';
mysyms = 'xosd';
tl = linspace(0,100);
for Bloadc = 1:numel(exd.Bsource)
  figure(1);
  ts = []; Bs = [];
  for expc = 1:numel(exd.t.all{Bloadc})
    t = exd.t.all{Bloadc}{expc};
    B = exd.B.all{Bloadc}{expc};
    p1 = plot(t,B,[mycols(Bloadc),mysyms(expc)]);
    p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    ts = [ts,t]; Bs = [Bs,B];
  end
  tsc{Bloadc} = ts;
  bsc{Bloadc} = Bs;
  my_model = @(a,tp) a(1)*exp(a(2)*tp);
  opts = statset('Display','iter','TolFun',1e-10,'MaxIter',1000);
  mdl = fitnlm(ts,Bs,my_model,[mean(Bs) -0.0001],'Options',opts)
  curve = fit(ts',log(Bs'),'poly1');
  [P,S] = polyfit(ts,log(Bs),1);
%   plot(tl,my_model(mdl.Coefficients.Estimate,tl),mycols(Bloadc));
%   plot(tl,curve.a*exp(curve.b*tl),mycols(Bloadc));
  plot(tl,exp(polyval(P,tl)),mycols(Bloadc));
  figure(2);
  plot(ts,abs(Bs-exp(polyval(P,ts))),...
    [mycols(Bloadc),mysyms(expc)]);
%   pause
end
figure(1);
set(gca,'YScale','log');
figure(2);
set(gca,'YScale','log');
figure(4);

close(figure(3)); figure(3);
for Bloadc = 1:numel(exd.Bsource)
  t = exd.t.ave{Bloadc};
  B = exd.B.ave{Bloadc}
  Bmin = exd.B.min{Bloadc}
  Bmax = exd.B.max{Bloadc}
  han = errorbar(t,B,B-Bmin,Bmax-B);
  set(han,'Color',mycols(Bloadc));
  hold on
end
for Bloadc = 1:numel(exd.Bsource)
  t = [exd.t.all{Bloadc}{:}];
  B = [exd.B.all{Bloadc}{:}];
  plot(t,B,[mycols(Bloadc),mysyms(expc)]);
end
set(gca,'YScale','log');
legend('1.28e8','2.48e8','5.05e8','1.94e?');

return

%  Calculate an exponential fit for all data
[stds{1:4}] = deal(1);

%  Various proposed algebraic expression models (no des).
%  Gaussian model for the data
mymodel = @(a,B,t) a(1)*B*exd(a(2)*(t-(a(3)*B+a(4))).^2);
a0 = [2,-0.001,-0.002,0];
%  Exponential decay model
mymodel = @(a,B,t) a(1)*B*exd(a(2)*t);
a0 = [2,-0.001,-0.002];

my_weighted_util = @(a,stds) ...
  sum((mymodel(a,Bsource(1),tsc{1})-bsc{1}).^2./stds{1}.^2)+...
  sum((mymodel(a,Bsource(2),tsc{2})-bsc{2}).^2./stds{2}.^2)+...
  sum((mymodel(a,Bsource(3),tsc{3})-bsc{3}).^2./stds{3}.^2)+...
  sum((mymodel(a,Bsource(4),tsc{4})-bsc{4}).^2./stds{4}.^2);
[a] = fminsearch(@(a) my_weighted_util(a,stds),[2,-0.001,-0.001,0]);
ts = [];
res = [];
jmax = 20;
for j = 1:jmax
  for Bloadc = 1:4
    resave = [];
    figure(1);
    p1(Bloadc) = plot(tl,mymodel(a,Bsource(Bloadc),tl),mycols(Bloadc));
    figure(4); if Bloadc == 1, hold off; else, hold on; end
    plot(tsc{Bloadc},log(resc{Bloadc}),[mycols(Bloadc),mysyms(expc)]);
    figure(2);
    resc{Bloadc} = abs(mymodel(a,Bsource(Bloadc),tsc{Bloadc})-bsc{Bloadc});
    p2(Bloadc) = plot(tsc{Bloadc},resc{Bloadc},[mycols(Bloadc),mysyms(expc)]);
    ts = [ts,tsc{Bloadc}];
    res = [res,resc{Bloadc}];
    tsunique = unique(tsc{Bloadc});
    for tuc = 1:numel(tsunique)
      inds = tsc{Bloadc} == tsunique(tuc);
      nsunique(tuc) = sum(inds);
      stds{Bloadc}(inds) = 1./nsunique(tuc);
      resave(tuc) = mean(resc{Bloadc}(tsunique(tuc) == tsc{Bloadc}));
    end
    myp = polyfit(tsunique,resave,1);
    stds{Bloadc} = polyval(myp,tsc{Bloadc}).*stds{Bloadc};
  end
  [a] = fminsearch(@(a) my_weighted_util(a,stds),[2,-0.001,-0.001,0]);
  if j < jmax
    delete(p1); delete(p2);
  end
end
figure(1);
legend('1.28e8','2.48e8','5.05e8','1.94e9');

close(figure(3)); figure(3); hold on
for Bloadc = 1:4
  plot(tsc{Bloadc},stds{Bloadc},[mycols(Bloadc),mysyms(expc)]);
end
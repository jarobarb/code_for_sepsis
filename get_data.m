function [exd] = get_data(outlier_matrix,varargin)

  %  This is stores all the experimental data (from the "Bacteria
  %  data*.xlsx")...including any outliers.
  %  All data (including outliers):
  %
  %  outlier_choice-allows us to pick outliers.  Currently, 0 corresponds
  %    to assuming there are no outliers, 1 picks off 1 outlier, and 2
  %    picks off the 1st outlier and another additional one
  %  td-Time delay before measurements start
  %  Bsource-initial bacterial levels
  td = 0;
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  exd.Bsource=[128 248 505 1940];

  %  exp.t.all{bacterial load index}{experiment index}
  exd.t.all{1}{1} = td+[24 48 72 96]; %Corresponds with Bsource1
  exd.B.all{1}{1} = [54.8 2.2 .74 .74];%averages of the concentrations
  exd.t.all{1}{2} = td+[24 48 72 96];
  exd.B.all{1}{2}= [6.6 4.6 3.2 1.2];
  exd.t.all{1}{3} = td+[24 48 96];
  exd.B.all{1}{3}= [21 160 .4];
  exd.t.all{1}{4} = td+[24 48 96];
  exd.B.all{1}{4} = [59 41.8 .4];

  exd.t.all{2}{1}  = td+[4 8 24 48 72]; %Corresponds with Bsource2
  exd.B.all{2}{1} = [150 144 190 8 2.4]; %averages of the concentrations

  exd.t.all{3}{1} = td+[0   2 2.5  3   4   6   8 18 28]; %Corresponds with Bsource3
  exd.B.all{3}{1} = [150 150 144 50 190 150 190 50 50]; %averages of the concentrations %outlier removed
  exd.t.all{3}{2} = td+[0 2 4];
  exd.B.all{3}{2}= [190 190 800];

  exd.t.all{4}{1} = td+[2 4 5 6 8 48]; %Corresponds with Bsource4
  exd.B.all{4}{1} = [800 190 50 190 800 140]; %averages of the concentrations %outlier removed
  exd.t.all{4}{2} = td+[2 4 6 48];
  exd.B.all{4}{2} = [150 150 190 150];
  exd.t.all{4}{3} = td+[2 4];
  exd.B.all{4}{3} = [190 110];
  exd.t.all{4}{4} = td+[2];
  exd.B.all{4}{4} = [110];

  %  Initialize "outlier" cells
  for Bloadc = 1:numel(exd.t.all)
    for expc = 1:numel(exd.t.all{Bloadc})
      exd.outlier{Bloadc}{expc} = false(size(exd.t.all{Bloadc}{expc}));
    end
  end

  %  Various choices for outliers
  for rc = 1:size(outlier_matrix,1)
    exd.outlier{outlier_matrix(rc,1)}{outlier_matrix(rc,2)}...
      (outlier_matrix(rc,3)) = true;
  end
  
  %  Now we go ahead and calculate average experimental values based on the
  %  above data (and outliers)
  exd.B.list = [];
  exd.t.list = [];
  exd.Bs.list = [];
  exd.B.lists = cell(1,numel(exd.Bsource));
  exd.t.lists = cell(1,numel(exd.Bsource));
  for Bloadc = 1:numel(exd.t.all)
    tunique = unique([exd.t.all{Bloadc}{:}]);
    exd.t.ave{Bloadc} = tunique;
    exd.B.sum{Bloadc} = 0*tunique;
    exd.B.lsum{Bloadc} = 0*tunique;
    exd.B.min{Bloadc} = 0*tunique; exd.B.min{Bloadc}(:) = inf;
    exd.B.max{Bloadc} = 0*tunique;
    exd.Bs.ave{Bloadc} = exd.B.sum{Bloadc}+exd.Bsource(Bloadc);
    exd.Bs.lave{Bloadc} = exd.B.lsum{Bloadc};
    exd.ns{Bloadc} = 0*tunique;
    exd.B.allo{Bloadc} = [];
    exd.Bs.allo{Bloadc} = [];
    exd.B.logallo{Bloadc} = [];
    exd.t.allo{Bloadc} = [];
    for expc = 1:numel(exd.t.all{Bloadc})
      for tc = 1:numel(exd.t.all{Bloadc}{expc})
        if ~exd.outlier{Bloadc}{expc}(tc)
          ind = find(exd.t.all{Bloadc}{expc}(tc) == tunique);
          exd.B.allo{Bloadc}(end+1) = exd.B.all{Bloadc}{expc}(tc);
          exd.Bs.allo{Bloadc}(end+1) = exd.Bsource(Bloadc);
          exd.t.allo{Bloadc}(end+1) = exd.t.all{Bloadc}{expc}(tc);
          exd.B.list(end+1) = exd.B.all{Bloadc}{expc}(tc);
          exd.t.list(end+1) = exd.t.all{Bloadc}{expc}(tc);
          exd.Bs.list(end+1) = exd.Bsource(Bloadc);
          exd.B.lists{expc}(end+1) = exd.B.all{Bloadc}{expc}(tc);
          exd.t.lists{expc}(end+1) = exd.t.all{Bloadc}{expc}(tc);
          exd.B.sum{Bloadc}(ind) = exd.B.sum{Bloadc}(ind)+...
            exd.B.all{Bloadc}{expc}(tc);
          exd.B.lsum{Bloadc}(ind) = exd.B.lsum{Bloadc}(ind)+...
            log(exd.B.all{Bloadc}{expc}(tc));
          exd.B.max{Bloadc}(ind) = max(exd.B.max{Bloadc}(ind),...
            exd.B.all{Bloadc}{expc}(tc));
          exd.B.min{Bloadc}(ind) = min(exd.B.min{Bloadc}(ind),...
            exd.B.all{Bloadc}{expc}(tc));
          exd.ns{Bloadc}(ind) = exd.ns{Bloadc}(ind)+1;
        end
      end
    end
    exd.B.ave{Bloadc} = exd.B.sum{Bloadc}./exd.ns{Bloadc};
    exd.B.logallo{Bloadc} = log(exd.B.allo{Bloadc});
    exd.B.lave{Bloadc} = exd.B.lsum{Bloadc}./exd.ns{Bloadc};
    exd.B.lave{Bloadc} = exp(exd.B.lave{Bloadc});
    if any(isnan(exd.B.ave{Bloadc}))
      inds = ~isnan(exd.B.ave{Bloadc});
      exd.B.ave{Bloadc} = exd.B.ave{Bloadc}(inds);
      exd.Bs.ave{Bloadc} = exd.Bs.ave{Bloadc}(inds);
      exd.t.ave{Bloadc} = exd.t.ave{Bloadc}(inds);
    end
  end
  
  %%  Remove any empty info
  %  For now, we just remove Bsource for any set of outliers that
  %  effectively removes all measurements for that Bsource.  Note that this
  %  is wired so that 1940 must be last in the list of experimental data.
  nBsource = numel(exd.Bsource);
  Bsourceinds = [];
  for Bsc = 1:nBsource
    if ~all(exd.ns{Bsc} == 0)
       Bsourceinds = [Bsourceinds,Bsc];
    end
  end
  exd.Bsource = exd.Bsource(Bsourceinds);
  exd.t.all = {exd.t.all{Bsourceinds}};
  exd.t.allo = {exd.t.allo{Bsourceinds}};
  exd.B.allo = {exd.B.allo{Bsourceinds}};
  exd.B.ave = {exd.B.ave{Bsourceinds}};
  exd.t.ave = {exd.t.ave{Bsourceinds}};
  exd.outlier_matrix = outlier_matrix;

end
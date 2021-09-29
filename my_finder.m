function [fvals,sol] = my_finder(varargin)

  fvals = varargin{2};
  for myc = 1:numel(fvals)
    [ystar,sol] = run_steady_state(varargin{1},varargin{2}(myc),...
      varargin{3},varargin{4}(myc),varargin{5:end});
    fvals(myc) = my_classifier(sol);%(ystar(1) > 1e-6) + (ystar(4) > 1);
    if nargout > 1
      for vac = 1:2:numel(varargin)
        varargin{vac} = strrep(varargin{vac},'Bsource_in','Bsource_ins');
      end
      [~,sol] = model_scaled_2018(varargin{1},varargin{2}(myc),...
        varargin{3},varargin{4}(myc),varargin{5:end});
    end
  end

end
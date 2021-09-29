function fval = my_classifier(sol)
  try
    septic_time = sol.extdata.varargin{3}.fint;
  catch
    septic_time = inf;
  end
  yinf = sol.y(:,end);
  fval = (yinf(1) > sol.extdata.varargin{3}.Binfasep) + ...
    (yinf(4) > sol.extdata.varargin{3}.dinfasep);
  %  High levels of bacteria, make it septic regardless of damage levels
  if (yinf(1) > sol.extdata.varargin{3}.Binfasep)
    fval = 2;
  end
  y = sol.y;  B = y(1,:); maxB = max(B);
  
  if sol.x(end) > septic_time
    yfin = deval(sol,septic_time);
    if yfin(1) > sol.extdata.varargin{3}.Bfinrsep*maxB, fval = 2; end
  end
%   fprintf('%30.20g\n',sol.extdata.varargin{3}.Bfinrsep);
%   end
end
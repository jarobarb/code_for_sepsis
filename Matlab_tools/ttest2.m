function [p,CIl,CIu] = ttest2(x,y,varargin)

  %  [p,CI] = ttest2(x,y)
  %  x and y must have the same number of columns (which could be 1).  The
  %    program gives the p-value associated with the probability that
  %    mean(x)-mean(y) does not equal zero.
  
  %  Default values for confidence interval stuff
  pCI = [0.05];
  [nx,nc] = size(x);
  [ny,nc] = size(y);
  mx = mean(x,1);
  my = mean(y,1);
  vx = var(x);
  vy = var(y);
	for vac = 1:2:numel(varargin)
		eval([varargin{vac},' = varargin{vac+1};']);
	end
  tstat = (mx-my)./sqrt(vx./nx+vy./ny);
  df = (vx./nx+vy./ny).^2./((vx./nx).^2./(nx-1)+(vy./ny).^2./(ny-1));
  
  CIl = []; CIu = [];
  
  for cc = 1:nc
    %  Two-tailed test/hypothesis is mean(x)-mean(y) is not equal to zero
    %  (vs just greater or just less than)
    p(cc) = 2*(1-tdf(abs(tstat(cc)),df(cc)));
  end
  
  for cc = 1:nc
    for CIc = 1:numel(pCI)
      dx = fzero(@(dx) pCI(CIc)-...
        2*(1-tdf(abs(dx./sqrt(vx(cc)/nx+vy(cc)/ny)),df(cc))),0);
      CIl(CIc,cc) = (mx(cc)-my(cc))-dx;
      CIu(CIc,cc) = (mx(cc)-my(cc))+dx;
    end
  end

end
function [yf,inds] = jared_lowess(x,y,varargin)

  %  default values
  drx = 0.05;
  num_repeats = 5;
  mad_fact = 6;
  %  Which outliers to throw away, positive (p), negative (n), or both (b)
  pos_or_neg_outliers = 'p';
  
  %  redefine default values
  for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  if length(who('dx')) == 0, dx = drx*(max(x)-min(x)); end
  yf = zeros(size(y));
  reg = ones(size(y));
  regws = reg;
  
  for j = 1:(num_repeats+1)
    
    not_adjusted = false(size(x));
    
    for xc = 1:length(x)
      
      inds = ((x(xc)-dx)<=x) & ((x(xc)+dx)>=x);
      tmpmax = max(abs(x(xc)-x(inds)));
      if tmpmax == 0
        ws = regws(inds);
      else
        ws = regws(inds).*(1-(abs(x(xc)-x(inds))./tmpmax).^3).^3;
      end
      if sum(ws) == 0
        warning(['No way to remove outlier using nearest neighbor(s).',...
          '  Using linear interpolation/extrapolation!']);
        not_adjusted(xc) = true;
%         if max(x(inds)) < max(x)
%           [~,indr] = min((x-x(xc)).*(x>(x(xc)+dx)).*(regws ~= 0)+...
%             (max(x)-min(x))*((x<=(x(xc)+dx))+(regws == 0)));
%         else
%           indr = [];
%         end
%         if min(x(inds)) > min(x)
%           [~,indl] = min((x(xc)-x).*(x<(x(xc)-dx)).*(regws ~= 0)+...
%             (max(x)-min(x))*((x>=(x(xc)-dx))+(regws == 0)));
%         else
%           indl = [];
%         end
%         
%         if ~isempty(indl) && ~isempty(indr)
%           yf(xc) = interp1(x([indl,indr]),y([indl,indr]),x(xc),'linear');
%         else
%           yf(xc) = mean(y([indr,indl]));
%         end
      else
        yf(xc) = sum(y(inds).*ws)./sum(ws);
      end
      reg(xc) = yf(xc)-y(xc);
      
    end
    
    %  Occasionally there are outliers with only fellow outliers as
    %  neighbors...define their values using linear
    %  interpolation/extrapolation
    yf(not_adjusted) = interp1(x(~not_adjusted),yf(~not_adjusted),...
      x(not_adjusted),'linear','extrap');
    
    if pos_or_neg_outliers == 'p'
      reg = -reg;
    elseif pos_or_neg_outliers == 'b';
      reg = abs(reg);
    end
       
    if j ~= (num_repeats+1)
      mad = median(abs(reg));
      regws = (1-(abs(reg)./mad./mad_fact).^2).^2.*(reg < (mad*mad_fact));
    else
      inds = (regws == 0);
    end
    
  end
    
  return;
  
  
  tmpx = xs(:);
  tmpx = [tmpx(1)*ones(n-1,1);tmpx;tmpx(end)*ones(n-1,1)];
  
  tmpy = ys(:);
  tmpy = [tmpy(1)*ones(n-1,1);tmpy;tmpy(end)*ones(n-1,1)];
  
  tmpxmat = repmat(tmpx,1,2*n-1);
  tmpymat = repmat(tmpy,1,2*n-1);
  
  for nc = 1:2*n-1
    tmpxmat(:,nc) = circshift(tmpx,nc-n);
    tmpymat(:,nc) = circshift(tmpy,nc-n);
    ds(:,nc) = abs(tmpxmat(:,nc)-tmpx);
  end
    
  maxd = max(ds,[],2);
  ws = (1-bsxfun(@times,ds.^3,1./maxd.^3)).^3;
  ws = ws(n:end-n+1,:);
  
  for nc = 1:n
    ws(nc,n+nc:end) = 0;
    ws(end-nc+1,1:n-nc) = 0;
  end
  
  yfilt = reshape(sum(tmpymat(n:end-n+1,:).*ws,2)./sum(ws,2),size(ys));
%   yreg = abs(ys-yfilt);
%   
%   ws2 = zeros(size(ws));
%   
%   for j = 1:5
%     
%     mad = median(abs(yreg));
%     wvec = (1-(yreg/6/mad).^2).^2.*(abs(yreg)<(6*mad));
%     for nc = 1:n
%       ws2(1:end-(n-nc),nc) = wvec((n-nc+1):end);
%       ws2((n-nc+1):end,nc) = wvec(1:end-(n-nc));
%     end
%     
%     yfilt = reshape(sum(tmpymat(n:end-n+1,:).*(ws.*ws2),2)./...
%       sum(ws.*ws2,2),size(ys));
%     yreg = ys-yfilt;
%     
%   end
  
end
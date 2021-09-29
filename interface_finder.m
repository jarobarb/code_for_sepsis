function [xg,yg,fg,fgo] = interface_finder(varargin)

  pr.ini_ref = 2;
  pr.nref = 2;
  pr.xl = 0;
  pr.xu = 1;
  pr.yl = 0;
  pr.yu = 1;
  pr.f = @(x,y) ((x.^2+y.^2) < 0.8^2)+((x-1).^2+(y-1).^2 < 0.8^2);
  pr.noutputs = 1;
  retest_boundary = true;
  test_all = false;
  pr.parcomp = true;
  
  for vac = 1:2:numel(varargin)
    pr.(varargin{vac}) = varargin{vac+1};
  end
  
  if test_all, ini_ref = (pr.ini_ref-1)*2^pr.nref+1;
  else, ini_ref = pr.ini_ref; end
  
  xvec = linspace(pr.xl,pr.xu,ini_ref);
  yvec = linspace(pr.yl,pr.yu,ini_ref);
  [xg,yg] = ndgrid(xvec,yvec);
  fg = xg;
  fgo = cell(size(fg));
  plot(xg,yg,'sr');
  
  if pr.parcomp
    p = gcp();
  end
  if ~pr.parcomp
    for myc = 1:numel(xg)
      [tmp{1:pr.noutputs}] = pr.f(xg(myc),yg(myc));
      fg(myc) = tmp{1};
      fgo{myc} = tmp{2:end};
    end
  else
    for myc = 1:numel(xg)
      my_jobs(myc) = parfeval(p,pr.f,pr.noutputs,xg(myc),yg(myc));
    end
    for myc = 1:numel(xg)
      [job_index,tmp{1:pr.noutputs}] = fetchNext(my_jobs);
      fg(job_index) = tmp{1};
      fgo{job_index} = tmp{2:end};
    end
  end
  if test_all, return; end
%   fg2 = pr.f(xg,yg);
  
  for rc = 1:pr.nref
    curr_ref = (pr.ini_ref-1)*2.^rc+1;
    xgprev = xg;
    ygprev = yg;
    fgprev = fg;
    fgoprev = fgo;
    
    xvec = linspace(pr.xl,pr.xu,curr_ref);
    yvec = linspace(pr.yl,pr.yu,curr_ref);
    
    [xg,yg] = ndgrid(xvec,yvec);
%     plot(xg,yg,'sr');
    fg(1:2:curr_ref,1:2:curr_ref) = fgprev;
    fgo = cell(size(fg));
    [fgo{1:2:curr_ref,1:2:curr_ref}] = deal(fgoprev{:});
    
    %  "Centered" midpoints
    for xt1c = 2:2:curr_ref-1
      for yt1c = 2:2:curr_ref-1
        if (fg(xt1c-1,yt1c-1) == fg(xt1c-1,yt1c+1)) && ...
            (fg(xt1c-1,yt1c-1) == fg(xt1c+1,yt1c-1)) && ...
            (fg(xt1c-1,yt1c-1) == fg(xt1c+1,yt1c+1))
          fg(xt1c,yt1c) = fg(xt1c-1,yt1c-1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,xt1c,yt1c,[1,1;1,-1;-1,1;-1,-1]);
          end
        else
          fg(xt1c,yt1c) = nan;%pr.f(xg(xt1c,yt1c),yg(xt1c,yt1c));
        end
      end
    end
    
    tmpind = find(isnan(fg))';
    tmpf = 0*tmpind;
    tmpfo = cell(size(tmpf));
    if pr.parcomp
      for linc = 1:numel(tmpind)
        ind = tmpind(linc);
        my_jobs(linc) = parfeval(p,pr.f,pr.noutputs,xg(ind),yg(ind));
      end
      for linc = 1:numel(tmpind)
        [job_index,tmp{1:pr.noutputs}] = fetchNext(my_jobs);
        tmpf(job_index) = tmp{1};
        tmpfo{job_index} = tmp{2:end};
      end
    else
      for linc = 1:numel(tmpind)
        ind = tmpind(linc);
        [tmp{1:pr.noutputs}] = pr.f(xg(ind),yg(ind));
        tmpf(linc) = tmp{1};
        tmpfo{linc} = tmp{2:end};
      end
    end
    fg(tmpind) = tmpf;
    [fgo{tmpind}] = deal(tmpfo{:});
    
    %  Boundaries of the domain
    for yt2xc = 2:2:curr_ref-1
      if (fg(1,yt2xc-1) == fg(1,yt2xc+1)) && ...
          (fg(1,yt2xc-1) == fg(2,yt2xc))
        if retest_boundary
          fg(1,yt2xc) = nan;%fg(1,yt2xc-1);
        else
          fg(1,yt2xc) = fg(1,yt2xc-1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,1,yt2xc,[0,1;0,-1;1,0]);
          end
        end
      else
        fg(1,yt2xc) = nan;%pr.f(xg(1,yt2xc),yg(1,yt2xc));
      end
      if (fg(end,yt2xc-1) == fg(end,yt2xc+1)) && ...
          (fg(end,yt2xc-1) == fg(end-1,yt2xc))
        if retest_boundary
          fg(end,yt2xc) = nan;%fg(end,yt2xc-1);
        else
          fg(end,yt2xc) = fg(end,yt2xc-1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,size(fg,1),yt2xc,[0,1;0,-1;-1,0]);
          end
        end
      else
        fg(end,yt2xc) = nan;%pr.f(xg(end,yt2xc),yg(end,yt2xc));
      end
    end
    for xt2yc = 2:2:curr_ref-1
      if (fg(xt2yc-1,1) == fg(xt2yc+1,1)) && ...
          (fg(xt2yc-1,1) == fg(xt2yc,2))
        if retest_boundary
          fg(xt2yc,1) = nan;%fg(xt2yc-1,1);
        else
          fg(xt2yc,1) = fg(xt2yc-1,1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,xt2yc,1,[-1,0;1,0;0,1]);
          end
        end
      else
        fg(xt2yc,1) = nan;%pr.f(xg(xt2yc,1),yg(xt2yc,1));
      end
      if (fg(xt2yc-1,end) == fg(xt2yc+1,end)) && ...
          (fg(xt2yc-1,end) == fg(xt2yc,end-1))
        if retest_boundary
          fg(xt2yc,end) = nan;%fg(xt2yc-1,end);
        else
          fg(xt2yc,end) = fg(xt2yc-1,end);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,xt2yc,size(fg,2),[-1,0;1,0;0,-1]);
          end
        end
      else
        fg(xt2yc,end) = nan;%pr.f(xg(xt2yc,end),yg(xt2yc,end));
      end
    end
    
    %  Midpoints at "x" edges
    for xt2xc = 3:2:curr_ref-2
      for yt2xc = 2:2:curr_ref-1
        if (fg(xt2xc,yt2xc-1) == fg(xt2xc,yt2xc+1)) && ...
            (fg(xt2xc,yt2xc-1) == fg(xt2xc-1,yt2xc)) && ...
            (fg(xt2xc,yt2xc-1) == fg(xt2xc+1,yt2xc))
          fg(xt2xc,yt2xc) = fg(xt2xc,yt2xc-1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,xt2xc,yt2xc,[0,1;0,-1;1,0;-1,0]);
          end
        else
          fg(xt2xc,yt2xc) = nan;%pr.f(xg(xt2xc,yt2xc),yg(xt2xc,yt2xc));
        end
      end
    end
    
    %  Midpoints at "y" edges
    for xt2yc = 2:2:curr_ref-1
      for yt2yc = 3:2:curr_ref-2
        if (fg(xt2yc,yt2yc-1) == fg(xt2yc,yt2yc+1)) && ...
            (fg(xt2yc,yt2yc-1) == fg(xt2yc-1,yt2yc)) && ...
            (fg(xt2yc,yt2yc-1) == fg(xt2yc+1,yt2yc))
          fg(xt2yc,yt2yc) = fg(xt2yc,yt2yc-1);
          if pr.noutputs > 1
            fgo = average_fgo(fgo,xt2yc,yt2yc,[0,1;0,-1;1,0;-1,0]);
          end
        else
          fg(xt2yc,yt2yc) = nan;%pr.f(xg(xt2yc,yt2yc),yg(xt2yc,yt2yc));
        end
      end
    end
    
    tmpind = find(isnan(fg))';
    tmpf = 0*tmpind;
    tmpfo = cell(size(tmpf));
    if pr.parcomp
      for linc = 1:numel(tmpind)
        ind = tmpind(linc);
        my_jobs(linc) = parfeval(p,pr.f,pr.noutputs,xg(ind),yg(ind));
      end
      for linc = 1:numel(tmpind)
        [job_index,tmp{1:pr.noutputs}] = fetchNext(my_jobs);
        tmpf(job_index) = tmp{1};
        tmpfo{job_index} = tmp{2:end};
      end
    else
      for linc = 1:numel(tmpind)
        ind = tmpind(linc);
        [tmp{1:pr.noutputs}] = pr.f(xg(ind),yg(ind));
        tmpf(linc) = tmp{1};
        tmpfo{linc} = tmp{2:end};
      end
    end
    fg(tmpind) = tmpf;
    [fgo{tmpind}] = deal(tmpfo{:});
        
%     contour(xg,yg,smooth_matrix(double(fg),0,0.5),[0.5,1.5],'b')
    indsn = isnan(fg);
    inds0 = (fg >= -0.5) & (fg < 0.5) & ~indsn;
    inds1 = (fg >= 0.5) & (fg < 1.5) & ~indsn;
    inds2 = (fg >= 1.5) & (fg < 2.5) & ~indsn;
    
    plot(xg(indsn),yg(indsn),'ko',xg(inds0),yg(inds0),'r<',...
      xg(inds1),yg(inds1),'b+',xg(inds2),yg(inds2),'g>');
%     surf(xg,yg,smooth_matrix(double(fg),0,0.5))
%     view(2)
%     colorbar
%     shading interp
%     ths = linspace(0,pi/2);
%     hold on
%     plot(0.8*cos(ths),0.8*sin(ths),'k');
    
    pause(0.01)
    
  end

end

function fgo = average_fgo(fgo,xi,yi,list)
  fn = fieldnames(fgo{xi+list(1,1),yi+list(1,2)});
  fn = setdiff(fn,'struc');
  for fc = 1:numel(fn)
    fgo{xi,yi}.(fn{fc}) = 0;
    for lc = 1:size(list,1)
      fgo{xi,yi}.(fn{fc}) = fgo{xi,yi}.(fn{fc})+...
        fgo{xi+list(lc,1),yi+list(lc,2)}.(fn{fc})(end);
    end
    fgo{xi,yi}.(fn{fc}) = fgo{xi,yi}.(fn{fc})/size(list,1);
  end
  fgo{xi,yi}.struc = fgo{xi+list(1,1),yi+list(1,2)}.struc;
end
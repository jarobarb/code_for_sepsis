function zi = tetralininterp(x,y,tri,xi,extrapval)

  %  When outside the given tetrahedral mesh, a value of extrapval is
  %  automatically assigned

  orig_siz = size(x);
  %  Dimension (2d or 3d)
  dim = orig_siz(2);
  % Now xi may = [X,Y] or [X,Y,Z] where X, Y, and Z are matrices rather
  matrix_siz(1) = size(xi,1);
  matrix_siz(2) = size(xi,2)/dim;
  xiold = xi;
  xi = [];
  for dc = 1:dim
    tmp = xiold(1:matrix_siz(1),(dc-1)*matrix_siz(2)+[1:matrix_siz(2)]);
    xi = [xi,tmp(:)];
  end
  % than column vectors.  We 
  % Find the nearest triangle (t)
  [t,p] = tsearchn(x,double(tri),xi);
  nan_inds = isnan(t);
  k = dsearchn(x,double(tri),xi(isnan(t),:));
  
  for nc = find(nan_inds)'
    if size(x,2) == 3
      [~,tmp_ind] = min((x(:,1)-xi(nc,1)).^2+(x(:,2)-xi(nc,2)).^2+...
        (x(:,3)-xi(nc,3)).^2);
    else
      [~,tmp_ind] = min((x(:,1)-xi(nc,1)).^2+(x(:,2)-xi(nc,2)).^2);
    end
    inds = find(sum(tri == tmp_ind,2));
    [t(nc),p(nc,:)] = tsearchn(x,double(tri(inds,:)),xi(nc,:));
    if isnan(t(nc))
      k(find(nc==find(nan_inds))) = tmp_ind;
    else
      t(nc) = inds(t(nc));
    end
  end
  nan_inds2 = isnan(t);

  m1 = size(xi,1);
  zi(size(xi,1),size(y,2)) = 0;% = nan(size(xi,1),size(y,2));

  for i = 1:m1
    if ~isnan(t(i))
      zi(i,:) = p(i,:)*y(tri(t(i),:),:);
    end
  end
  if isstr(extrapval)
    zi(nan_inds,:) = y(k,:);
  else
    zi(nan_inds,:) = extrapval;
  end
  
  %  Reshape if need be
  if matrix_siz(2) > 1    
    zi = reshape(zi,matrix_siz(1),matrix_siz(2));
  end

end
function [graphs,vertices,y,iter,fail] = par_simplex(vertices,fhereiam,varargin)

if ~isempty(varargin)
  tol = varargin{1};
else
  tol = 1e-3;
end

done = 0;
maxits = 500;
if min(size(vertices)) == 1
  vertices = vertices(:)';
  vertices = [vertices;
    vertices.*(ones(numel(vertices),1)+100*tol*eye(numel(vertices)))];
end
ndim = min(size(vertices));

parcomp = true;
if parcomp
  p = gcp();

  for n = 1:ndim+1
    my_jobs(n) = parfeval(p,fhereiam,1,vertices(n,:)');
  end
  y = zeros(1,ndim+1);
  for n = 1:ndim+1
    [job_index,value] = fetchNext(my_jobs);
    y(job_index) = value;
    fprintf('Job #%d/%d done\n',job_index,ndim+1);
  end
else
  for n = 1:ndim+1
    y(n) = feval(fhereiam,vertices(n,:)');
    fprintf('%g %d/%d done.\n',y(n));
  end
end

psum = sum(vertices,1);
iter = 0;
graphs = zeros(2,0);
fail = 0;

while ~done
    fprintf('%g\n',min(y));
    orderedy = sort(y);
    [tf,ihi] = ismember(orderedy(ndim+1),y);
    [tf,inhi] = ismember(orderedy(ndim),y);
    [tf,ilo] = ismember(orderedy(1),y);
    rtol = 0;
    for ii=1:ndim+1
        for j=ii+1:ndim+1
            rtol = max(norm(vertices(ii,:)-vertices(j,:)),rtol);
        end
    end
%     rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)));
    if rtol < tol
        swap=y(1);
        y(1)=y(ilo);
        y(ilo)=swap;
        for n=1:ndim
            swap=vertices(1,n);
            vertices(1,n)=vertices(ilo,n);
            vertices(ilo,n)=swap;
        end
        return
    end
    if (iter >= maxits)
        fail = 1;
        return
    end
    iter = iter+2;
    [vertices,ytry,psum,y] = flipper(vertices,y,psum,ndim,fhereiam,ihi,-1);
    if ytry <= y(ilo)
        [vertices,ytry,psum,y] = flipper(vertices,y,psum,ndim,fhereiam,ihi,2);
    elseif ytry >= y(inhi)
        ysave = y(ihi);
        [vertices,ytry,psum,y] = flipper(vertices,y,psum,ndim,fhereiam,ihi,.5);
        if ytry >= ysave
            for ii = 1:ndim+1
                if ii ~= ilo
                    for j = 1:ndim
                        psum(j) = .5*(vertices(ii,j)+vertices(ilo,j));
                        vertices(ii,j) = psum(j);
                    end
                    my_jobs(ii) = parfeval(p,fhereiam,1,psum');                    
                end
            end
            for ii = 1:ndim+1
              if ii ~= ilo
                [job_index,value] = fetchNext(my_jobs);
                y(job_index) = value;
                fprintf('Job #%d/%d done\n',job_index,ndim+1);
              end
            end
            iter = iter+ndim;
            psum = sum(vertices,1);
        end
    else
        iter = iter - 1;
    end
    graphs = [graphs,[iter;min(y)]];
end

end

function [vertices,ytry,psum,y] = flipper(vertices,y,psum,ndim,f,ihi,fac)

  fac1 = (1.-fac)/ndim;
  fac2 = fac1-fac;
  for j = 1:ndim
      ptry(j) = psum(j)*fac1-vertices(ihi,j)*fac2;
  end
  ytry = feval(f,ptry');
  if ytry < y(ihi)
      y(ihi) = ytry;
      for j=1:ndim
          psum(j) = psum(j)-vertices(ihi,j)+ptry(j);
          vertices(ihi,j) = ptry(j);
      end
  end

end
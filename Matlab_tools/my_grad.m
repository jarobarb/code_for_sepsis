function [fdx,fdy,fdz] = my_grad(f,xyz,origtetra)

  dim = 2+isfield(xyz,'z');
  fdx = f; fdy = f; fdz = f;
  
  ndof = size(origtetra,2);
  Mx = sparse(ndof,ndof);
  My = sparse(ndof,ndof);
  Mz = sparse(ndof,ndof);
  
  switch dim
    case 2
      error('2d not done yet');
    case 3
      switch size(origtetra,2)
        %  Quadratic elements
        case 10          
          for tc = 1:size(origtetra,1)
            inds = origtetra(tc,:)';
            locx = xyz.x(inds); locy = xyz.y(inds); locz = xyz.z(inds);
            M = [locx.^2,locy.^2,locz.^2,...
              locx.*locy,locx.*locz,...
              locy.*locz,...
              locx,locy,locz,ones(size(inds))];
            Mx(:,1) = 2.*locx;
            Mx(:,4) = locy;
            Mx(:,5) = locz;
            Mx(:,7) = 1;
            My(:,2) = 2.*locy;
            My(:,4) = locx;
            My(:,6) = locz;
            My(:,8) = 1;
            Mz(:,3) = 2.*locz;
            Mz(:,5) = locx;
            Mz(:,6) = locy;
            Mz(:,9) = 1;
            
            as = M\f(inds);
            fdx(inds) = Mx*as;
            fdy(inds) = My*as;
            fdz(inds) = Mz*as;            
          end
          
      end
  end

end
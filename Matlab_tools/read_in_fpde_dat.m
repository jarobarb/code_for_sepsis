function [xyzs,rest,DT,tetra,tri,t0,origtetra,node_info] =...
	read_in_fpde_dat(varargin)

  %  Typical usage:
  %  [xyzs,rest,DT,tetra,tri,t0,origtetra,node_info] =...
	%    read_in_fpde_dat('filename','Stress_2d.dat');
	%  xyzs-[xs,ys,zs] coordinates for nodes for mesh
	%  rest-values at the nodes for solutions or any other field exported
	%    using the transfer mesh function (not complicated exports
	%    though...e.g. trying transfer(dx(u)) may or may not work
	%  DT-delaunay triangulation which may or may not be empty depending on
	%    do_delaunay option
	%  tetra-multiple rows, each row corresponding to a tetrahedron that
	%    makes up the mesh and the values in the rows being the index for the
	%    nodes that make up the tetrahedron (4 nodes)
	%  tri-multiple rows, each row corresponding to triangles on the
	%    tetrahedron that make up the mesh, see tetra
	%  t0-time (if appropriate)
	%  origtetra-the original elements which may have 4 nodes apiece (linear
	%    tetrahedral elements) 10 nodes apiece (quadratic elements) and some
	%    27 nodes apiece (cubic elements) meshes may even work
	%  node_info-has info regarding node classification like whether or not
	%    the node is in the interior, on one side of the material, at a
	%    corner, and perhaps some other info...used in tetra_exterior_plot.m
	

  %  Default filename for debugging
	%  3d example
  filename = 'Stress_2d.dat';
	%  2d example
	filename = 'Stress_2d_wedge_plane_stress_flat.dat';
  do_delaunay = false;
  ele_vertices_only = false;
  
  for vac = 1:2:length(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  fid = fopen(filename);
  
  %  Get time (if available) and number of dimensions
  str = fgetl(fid);
  t0 = [];
  while isempty(strfind(str,'DIMENSIONS'))
    if ~isempty(strfind(str,'Time'))
      t0 = sscanf(str(6:end),'%g');
    end
    str = fgetl(fid);    
  end
  dim = sscanf(str,'%g');
  
  %  Get order of elements
  while isempty(strfind(str,'DOF'))
    str = fgetl(fid);
  end
  ndof = sscanf(str,'%g');
  if dim == 3
    if ndof == 10
      order = 2;
    elseif ndof == 20
      order = 3;
    else
      order = 1;
    end
  else
    if ndof == 6
      order = 2;
    elseif (ndof == 10) || (ndof == 9)
      order = 3;
    else
      order = 1;
    end
  end
  
  %  Get number of variables and variable names
  while isempty(strfind(str,'VARIABLES'))
    str = fgetl(fid);    
  end
  nvar = sscanf(str,'%g');
  
  vc = 0;
  
  while vc < nvar
    str = fgetl(fid);
    [~,k] = ismember(str,'=');
    k = find(k);
    if ~isempty(k)
      rest.(strrep(str(1:k-2),' ','')) = str2num(str(k+1:end));
      nvar = nvar-1;
    else
      vc = vc+1;
      var_name{vc} = sscanf(str,'%s');
      var_name{vc} = strrep(strrep(strrep(var_name{vc},'(',''),')',''),...
        '"','');
      [~,k] = ismember(var_name{vc},'-');
      k = find(k);
      if ~isempty(k)
        var_name{vc} = var_name{vc}(1:k-2);
      end
      leading_num = str2num(var_name{vc}(1));
      my_num2str = 'abcdefghi';
      if ~isempty(leading_num)
        var_name{vc}(1) = my_num2str(leading_num);
      end
    end
  end
  
  %  Get number of regions
  while isempty(strfind(str,'REGIONS'))
    str = fgetl(fid);
  end
  nreg = sscanf(str,'%g');
  for rc = 1:nreg
    str = fgetl(fid);
    my_mat(rc,:) = sscanf(str,'%d %d %d %d %d');
    k = strfind(str,'"');
    my_reg{rc} = str(k(1)+1:k(2)-1);
  end
  tmp_reg_inds = my_mat(:,end-1);
  [reg_inds,tmp_inds] = unique(tmp_reg_inds);
  regs = {my_reg{tmp_inds}};
  nreg = max(reg_inds);
  
  for rc = 1:nreg
    
    %  Get nodes and info at/on nodes
    while isempty(strfind(str,'NODES'))
      str = fgetl(fid);
    end
    nnod = sscanf(str,'%g');
    
    my_scan_str = '';
    for ssc = 1:dim+nvar
      my_scan_str = [my_scan_str,'%f64 '];
    end
    my_scan_str = [my_scan_str,' : %d %d %d %d'];
    
    C = textscan(fid,my_scan_str,nnod);
    
    n = numel(C);
    node_info = [C{n-3:n}];
    xyzs.x = C{1};
    xyzs.y = C{2};
    if dim == 3, xyzs.z = C{3}; end
    for vc = 1:nvar
      rest.(var_name{vc}) = C{vc+dim};
    end
    rest.reg_name = regs{rc};
    if isfield(rest,'sxy')
%       if isfield(rest,'txy')
%         rest.sxy = rest.txy;
%         rest.syz = rest.tyz;
%         rest.sxz = rest.tzx;
%       end
      rest.smax(numel(rest.sx),1) = 0;
      rest.smid(numel(rest.sx),1) = 0;
      rest.smin(numel(rest.sx),1) = 0;
      rest.nmaxx(numel(rest.sx),1) = 0;
      rest.nmaxy(numel(rest.sx),1) = 0;
      rest.nmaxz(numel(rest.sx),1) = 0;
      rest.nmidx(numel(rest.sx),1) = 0;
      rest.nmidy(numel(rest.sx),1) = 0;
      rest.nmidz(numel(rest.sx),1) = 0;
      rest.nminx(numel(rest.sx),1) = 0;
      rest.nminy(numel(rest.sx),1) = 0;
      rest.nminz(numel(rest.sx),1) = 0;
      for n = 1:length(rest.sx)
        M = [rest.sx(n),rest.sxy(n),rest.sxz(n);...
          rest.sxy(n),rest.sy(n),rest.syz(n);...
          rest.sxz(n),rest.syz(n),rest.sz(n)];
        [V,D] = eig(M);
        rest.smax(n) = D(3,3);
        rest.nmaxx(n) = V(1,3);
        rest.nmaxy(n) = V(2,3);
        rest.nmaxz(n) = V(3,3);
        rest.smid(n) = D(2,2);
        rest.nmidx(n) = V(1,2);
        rest.nmidy(n) = V(2,2);
        rest.nmidz(n) = V(3,2);
        rest.smin(n) = D(1,1);
        rest.nminx(n) = V(1,1);
        rest.nminy(n) = V(2,1);
        rest.nminz(n) = V(3,1);
        if (mod(n,10000) == 0), fprintf('%d/%d\n',n,length(rest.sx)); end
      end
      rest.scrit = rest.Von_mises./rest.relvmmed;
      Tps = [22,30,50,70,78];
      Comp_strengths = [44,43,38,29,24];
      x = fminsearch(@(x) norm(Comp_strengths-...
        (x(1)-exp(x(2)*(Tps-x(3))))),[45,1,20]);
      scrit_temp = 1e6*(x(1)-exp(x(2)*(rest.Tp-x(3))));
      %     rest.scrit2 = (rest.Tp < 40).*scrit_temp+(rest.Tp > 60).*rest.scrit+...
      %       ((rest.Tp >= 40) & (rest.Tp <= 60)).*...
      %       (pchip([20,40,60,80],[0,0,1,1],rest.Tp).*rest.scrit+...
      %       (1-pchip([20,40,60,80],[0,0,1,1],rest.Tp)).*scrit_temp);
      rest.scrit2 = (rest.Tp < 40).*scrit_temp+(rest.Tp > 60).*rest.scrit+...
        ((rest.Tp >= 40) & (rest.Tp <= 60)).*...
        (min(rest.scrit,scrit_temp));
      %     rest.scrit2 = min(rest.scrit,scrit_temp);
    end
        
    str = fgetl(fid); str = fgetl(fid);
    nquadele = sscanf(str,'%g');
    npost_colon = 4+dim;
    
    my_scan_str = '%d ';
    for ssc = 1:ndof
      my_scan_str = [my_scan_str,'%d '];
    end
    my_scan_str = [my_scan_str,' :'];
    for ssc = 1:npost_colon
      my_scan_str = [my_scan_str,' %d'];
    end
    C = textscan(fid,my_scan_str,nquadele);
    
    origtetra = [C{2:1+ndof}];
        
    if isfield(rest,'AxR') & ~isfield(rest,'JxR')
      if ~isfield(rest,'VR')
        rest.VR = 0*rest.AxR;
        rest.VI = 0*rest.AxR;
      end
      rest.dxAxR = zeros(size(xyzs.x));
      rest.dxAyR = zeros(size(xyzs.x));
      rest.dyAxR = zeros(size(xyzs.x));
      rest.dyAyR = zeros(size(xyzs.x));
      rest.dxAxI = zeros(size(xyzs.x));
      rest.dxAyI = zeros(size(xyzs.x));
      rest.dyAxI = zeros(size(xyzs.x));
      rest.dyAyI = zeros(size(xyzs.x));
      rest.dxVR = zeros(size(xyzs.x));
      rest.dyVR = zeros(size(xyzs.x));
      rest.dxVI = zeros(size(xyzs.x));
      rest.dyVI = zeros(size(xyzs.x));
      my_count = zeros(size(xyzs.x));
      Mx = sparse(ndof,ndof);
      My = sparse(ndof,ndof);
      Mz = sparse(ndof,ndof);
      tic;
      fprintf('Calculating mag/scalar potential derivatives in region %s...',...
        num2str(rc));
      for ac = 1:size(origtetra,1)
        if dim == 3
          inds = origtetra(ac,:)';
          locx = xyzs.x(inds); locy = xyzs.y(inds); locz = xyzs.z(inds);
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
          if order == 3
            M = [M,...
              locx.^3,locy.^3,locz.^3,...
              locx.^2.*locy,locx.^2.*locz,...
              locy.^2.*locx,locy.^2.*locz,...
              locz.^2.*locx,locz.^2.*locy,...
              locx.*locy.*locz];
            Mx(:,11) = 3*locx.^2;
            My(:,12) = 3*locy.^2;
            Mz(:,13) = 3*locz.^2;
            Mx(:,14) = 2*locx.*locy;
            My(:,14) = locx.^2;
            Mx(:,15) = 2*locx.*locz;
            Mz(:,15) = locx.^2;
            Mx(:,16) = locy.^2;
            My(:,16) = 2*locx.*locy;
            My(:,17) = 2*locy.*locz;
            Mz(:,17) = locy.^2;
            Mx(:,18) = locz.^2;
            Mz(:,18) = 2*locx.*locz;
            My(:,19) = locy.^2;
            Mz(:,19) = 2*locy.*locz;
            Mx(:,20) = locy.*locz;
            My(:,20) = locx.*locz;
            Mz(:,20) = locx.*locy;
          end
          
          as = M\rest.AxR(inds);        
          rest.dxAxR(inds) = Mx*as;
          rest.dyAxR(inds) = My*as;
          as = M\rest.AyR(inds);
          rest.dxAyR(inds) = Mx*as;
          rest.dyAyR(inds) = My*as;

          as = M\rest.AxR(inds);        
          rest.dxAxI(inds) = Mx*as;
          rest.dyAxI(inds) = My*as;
          as = M\rest.AyR(inds);
          rest.dxAyI(inds) = Mx*as;
          rest.dyAyI(inds) = My*as;

          as = M\rest.VR(inds);
          rest.dxVR(inds) = Mx*as;
          rest.dyVR(inds) = My*as;
          as = M\rest.VI(inds);
          rest.dxVI(inds) = Mx*as;
          rest.dyVI(inds) = My*as;
          
          my_count(inds) = my_count(inds)+1;
        end
      end
      rest.dxAxR = rest.dxAxR./my_count;
      rest.dxAyR = rest.dxAyR./my_count;
      rest.dyAxR = rest.dyAxR./my_count;
      rest.dyAyR = rest.dyAyR./my_count;
      rest.dxAxI = rest.dxAxI./my_count;
      rest.dxAyI = rest.dxAyI./my_count;
      rest.dyAxI = rest.dyAxI./my_count;
      rest.dyAyI = rest.dyAyI./my_count;
      rest.dxVR = rest.dxVR./my_count;
      rest.dyVR = rest.dyVR./my_count;
      rest.dxVI = rest.dxVI./my_count;
      rest.dyVI = rest.dyVI./my_count;
      
      rest.BzR = rest.dxAyR-rest.dyAxR;
      rest.BzI = rest.dxAyI-rest.dyAxI;
      freq = 100;
      fprintf('Frequency = %g Hz.\n',freq);
      omega = freq*2*pi;
      rest.ExR = omega*rest.AxI-rest.dxVR;
      rest.EyR = omega*rest.AyI-rest.dyVR;
      rest.ExI = -omega*rest.AxR-rest.dxVI;
      rest.EyI = -omega*rest.AyR-rest.dyVI;
      fprintf('done in %g seconds.\n',toc);
    end
    
    if dim == 2
      
      if order == 1
        origeletri = [1,2,3];
      elseif order == 2
        origeletri = [1,6,5;...
          6,2,4;...
          6,4,5;...
          5,4,3];
      else
        my_hexagon = [4,5,6,7,8,9];
        tmpxs = mean(xyzs.x(origtetra(:,my_hexagon)),2);
        xyzs.x = [xyzs.x;tmpxs];
        tmpys = mean(xyzs.y(origtetra(:,my_hexagon)),2);
        xyzs.y = [xyzs.y;tmpys];
        fn = fieldnames(rest);
        for fnc = 1:numel(fn)
          if numel(rest.(fn{fnc})) == nnod
            tmprest = mean(rest.(fn{fnc})(origtetra(:,my_hexagon)),2);
            rest.(fn{fnc}) = [rest.(fn{fnc});tmprest];
          end
        end
        origtetra = [origtetra,nnod+[1:nquadele]'];
        node_info = [node_info;zeros(nquadele,4)];
        nnod = numel(xyzs.x);
        nquadele = size(origtetra,1);
%         meanxs = mean(xyzs.x(origtetra),2);
%         meanys = mean(xyzs.y(origtetra),2);
%         dist2center = (xyzs.x(origtetra)-meanxs).^2+...
%           (xyzs.y(origtetra)-meanys).^2;
%         [tmp,ind10] = min(dist2center,[],2);
%         lind10 = sub2ind(size(origtetra),[1:numel(ind10)]',ind10);
%         
%         [tmp,ind1] = max(dist2center,[],2);
%         lind1 = sub2ind(size(origtetra),[1:numel(ind1)]',ind1);
%         dist21 = (xyzs.x(origtetra)-xyzs.x(origtetra(lind1))).^2+...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind1))).^2;
%         [tmp,ind2] = max(dist21,[],2);
%         lind2 = sub2ind(size(origtetra),[1:numel(ind2)]',ind2);
%         dist212 = abs((xyzs.x(origtetra)-xyzs.x(origtetra(lind1))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind2)))-...
%           (xyzs.x(origtetra)-xyzs.x(origtetra(lind2))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind1))));
%         [tmp,ind3] = max(dist212,[],2);
%         lind3 = sub2ind(size(origtetra),[1:numel(ind3)]',ind3);
%         
%         [tmp,ind69] = mink(dist212,4,2);
%         onenumv = [1:numel(ind1)]';
%         onetenm = [1:10].*ones(numel(ind1),1);
%         logind69 = (dist212 <= max(tmp,[],2)) & (onetenm ~= ind1) & ...
%           (onetenm ~= ind2);
%         [a,b] = ind2sub(size(origtetra'),find(logind69'));
%         ind69 = reshape(a,2,numel(ind1))';
%         linds69 = sub2ind(size(dist21),[onenumv,onenumv],ind69);
%         [tmp,ind6] = min(dist21(ind69),[],2);
%         lind6 = sub2ind(size(ind69),onenumv,ind6);
%         ind9 = 3-ind6;
%         lind9 = sub2ind(size(ind69),onenumv,ind9);
%         ind6 = ind69(lind6);
%         ind9 = ind69(lind9);
%         lind6 = sub2ind(size(origtetra),[1:numel(ind6)]',ind6);
%         lind9 = sub2ind(size(origtetra),[1:numel(ind9)]',ind9);
%         
%         dist22 = (xyzs.x(origtetra)-xyzs.x(origtetra(lind2))).^2+...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind2))).^2;
%         dist223 = abs((xyzs.x(origtetra)-xyzs.x(origtetra(lind2))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind3)))-...
%           (xyzs.x(origtetra)-xyzs.x(origtetra(lind3))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind2))));
%         [tmp,ind47] = mink(dist223,4,2);
%         logind47 = (dist223 <= max(tmp,[],2)) & (onetenm ~= ind2) & ...
%           (onetenm ~= ind3);
%         [a,b] = ind2sub(size(origtetra'),find(logind47'));
%         ind47 = reshape(a,2,numel(ind2))';
%         [tmp,ind4] = min(dist22(ind47),[],2);
%         lind4 = sub2ind(size(ind47),onenumv,ind4);
%         ind7 = 3-ind4;
%         lind7 = sub2ind(size(ind47),onenumv,ind7);
%         ind4 = ind47(lind4);
%         ind7 = ind47(lind7);
%         lind4 = sub2ind(size(origtetra),[1:numel(ind4)]',ind4);
%         lind7 = sub2ind(size(origtetra),[1:numel(ind7)]',ind7);
%         
%         dist23 = (xyzs.x(origtetra)-xyzs.x(origtetra(lind2))).^2+...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind2))).^2;
%         dist231 = abs((xyzs.x(origtetra)-xyzs.x(origtetra(lind3))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind1)))-...
%           (xyzs.x(origtetra)-xyzs.x(origtetra(lind1))).*...
%           (xyzs.y(origtetra)-xyzs.y(origtetra(lind3))));
%         [tmp,ind58] = mink(dist231,4,2);
%         logind58 = (dist231 <= max(tmp,[],2)) & (onetenm ~= ind3) & ...
%           (onetenm ~= ind1);
%         [a,b] = ind2sub(size(origtetra'),find(logind58'));
%         ind58 = reshape(a,2,numel(ind3))';
%         [tmp,ind5] = min(dist23(ind58),[],2);
%         lind5 = sub2ind(size(ind58),onenumv,ind5);
%         ind8 = 3-ind5;
%         lind8 = sub2ind(size(ind58),onenumv,ind8);
%         ind5 = ind58(lind5);
%         ind8 = ind58(lind8);
%         lind5 = sub2ind(size(origtetra),[1:numel(ind5)]',ind5);
%         lind8 = sub2ind(size(origtetra),[1:numel(ind8)]',ind8);
%         
%         oldtetra = origtetra;
%         origtetra = origtetra([lind1,lind2,lind3,lind4,lind5,lind6,lind7,lind8,lind9,lind10]);
        
        origeletri = [1,8,6;...
          5,3,7;...
          9,4,2;...
          8,6,10;...
          8,5,10;...
          5,7,10;...
          7,4,10;...
          4,9,10;...
          9,6,10];
%         origeletri = [1,2,3];
      end        
      
      tetra = reshape(origtetra(:,origeletri')',3,[])';      
      tri = tetra;
      
    else
      %  Assume maximum 10 mini tetrahedra for each quad tetrahedra
      %     max_tetra_per_quad_tetra = 26;
      %     tetra = -ones(size(quadtetra,1)*max_tetra_per_quad_tetra,4);
      %     for j = 1:size(quadtetra,1)
      %       quadeletri = delaunay(xyzs.x(quadtetra(j,:)),xyzs.y(quadtetra(j,:)),...
      %         xyzs.z(quadtetra(j,:)));
      %       if size(quadeletri,1) > max_tetra_per_quad_tetra
      %         error('Try increasing max_tetra_per_quad_tetra.');
      %       end
      %       tetra((j-1)*max_tetra_per_quad_tetra+[1:size(quadeletri,1)],:) =...
      %         reshape(quadtetra(j,quadeletri),[],4);
      %       if mod(j,10000) == 0, fprintf('%g/%g\n',j,size(quadtetra,1)); end
      %     end
      %     tetra = tetra(tetra(:,1)>0,:);
      %     quadeletri = [10,6,5,2;...
      %       4,9,8,10;...
      %       3,5,6,10;...
      %       7,8,1,6;...
      %       10,6,8,2;...
      %       2,8,7,6;...
      %       9,8,2,10];
      if ele_vertices_only
        quadeletri = [1,2,3,4];
      else
        quadeletri = [2,7,9,5;...
          6,10,3,5;...
          9,8,10,4;...
          1,6,8,7;...
          7,5,9,8;...
          5,9,10,8;
          5,7,6,8;
          5,10,6,8];
      end
      
      if ndof == 20
        if ele_vertices_only
          origeletri = [1,2,3,4];
        else
          tri1 = [14,15,16,4,17,18,19,8,9,10];
          tri2 = [1,13,6,8,20,12,7,14,19,18];
          tri3 = [12,5,3,10,11,6,20,18,17,16];
          tri4 = [7,2,11,9,5,20,13,19,15,17];
          tmptri = [tri1(quadeletri);tri2(quadeletri);...
            tri3(quadeletri);tri4(quadeletri)];
          ntris = size(tmptri,1);
          my_include = true(ntris);
          for c1 = 1:ntris
            for c2 = c1+1:ntris
              my_include(c1,c2) = numel(union(tmptri(c1,:),...
                tmptri(c2,:)))>4;
            end
          end
          origeletri = tmptri(all(my_include,2),:);
        end
      elseif ndof == 4
        origeletri = [1,2,3,4];
      else
        origeletri = quadeletri;
      end
      for oetc = 1:size(origeletri,1)
        rear = delaunay([xyzs.x(origtetra(1,origeletri(oetc,:))),...
          xyzs.y(origtetra(1,origeletri(oetc,:))),...
          xyzs.z(origtetra(1,origeletri(oetc,:)))]);
        origeletri(oetc,:) = origeletri(oetc,rear);
      end
      
      tetra = reshape(origtetra(:,origeletri')',4,[])';
      %     for tetrac = 1:size(tetra,1)
      %       tetra(tetrac,:) = tetra(tetrac,delaunay(xyzs.x(tetra(tetrac,:)),...
      %         xyzs.y(tetra(tetrac,:)),xyzs.z(tetra(tetrac,:))));
      %       if (mod(tetrac,1000)) == 0
      %         fprintf('%g/%g\n',tetrac,size(tetra,2));
      %       end
      %     end
      %     tetra = reshape(quadtetra(:,quadeletri'),[],4);
      tetratotri = [1,2,3;1,4,3;1,2,4;2,3,4];
      if ele_vertices_only
        tri = reshape(tetra(:,tetratotri')',3,[])';
      else
%         max_tri_per_orig_ele = 2*25;
%         tri = -ones(size(origtetra,1)*max_tri_per_orig_ele,3);
%         for j = 1:size(origtetra,1)
%           tmptri = convhull(xyzs.x(origtetra(j,:)),xyzs.y(origtetra(j,:)),...
%             xyzs.z(origtetra(j,:)));
%           if size(tmptri,1) > max_tri_per_orig_ele
%             error('Try increasing max_tri_per_quad');
%           end
%           tri((j-1)*max_tri_per_orig_ele+[1:size(tmptri,1)],:) = ...
%             reshape(origtetra(j,tmptri),[],3);
%         end
%         tri = tri(tri(:,1)>0,:);
        tri = reshape(tetra(:,tetratotri')',3,[])';
      end
    end
    
    if do_delaunay
      if dim == 2
        DT = DelaunayTri(xyzs.x,xyzs.y);
      else
        DT = DelaunayTri(xyzs.x,xyzs.y,xyzs.z);
      end
    else
      DT = {};
    end
    
    if nreg >= 1
      xyzss{rc} = xyzs;
      rests{rc} = rest;
      DTs{rc} = DT;
      tetras{rc} = tetra;
      tris{rc} = tri;
      t0s{rc} = t0;
      origtetras{rc} = origtetra;
      node_infos{rc} = node_info;
    end
    
  end
  
  if nreg > 1
    xyzs = xyzss;
    rest = rests;
    DT = DTs;
    tetra = tetras;
    tri = tris;
    t0 = t0s;
    origtetra = origtetras;
    node_info = node_infos;
  else
    xyzs = xyzss{1};
    rest = rests{1};
    DT = DTs{1};
    tetra = tetras{1};
    tri = tris{1};
    t0 = t0s{1};
    origtetra = origtetras{1};
    node_info = node_infos{1};
  end
    
  fclose(fid);

end
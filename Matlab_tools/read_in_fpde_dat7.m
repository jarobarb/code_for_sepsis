function [xyz,rest,origtetra,tetra,ni] =...
	read_in_fpde_dat7(varargin)

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
  while ~contains(str,'DIMENSIONS')
    if contains(str,'Time')
      t0 = sscanf(str(6:end),'%g');
    end
    str = fgetl(fid);    
  end
  dim = sscanf(str,'%g');
  coords = 'xyz'; coords = coords(1:dim);
  
  %  Get order of elements
  while ~contains(str,'CELLDOF')
    str = fgetl(fid);
  end
  dof = sscanf(str,'%g');
  if dim == 3
    if dof == 10
      order = 2;
    elseif dof == 20
      order = 3;
    else
      order = 1;
    end
  else
    if dof == 6
      order = 2;
    elseif (dof == 10) || (dof == 9)
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
  
  for rc = 1:1%nreg
    
    %  Get cell corner nodes
    while isempty(strfind(str,'NODES'))
      str = fgetl(fid);
    end
    ncorn = sscanf(str,'%g');
    
    %  Corner/joint node index, x coord, y coord, the node type, 
    %  (0=interior; 1=joint; 2=edge; 3=face; 4=exterior), boundary
    %  identifier (region number, joint number, edge number or face number)
    %  coefficient index (i.e. global node index...useful)
    if dim == 2
      C = textscan(fid,'%d %f64 %f64 %d %d %d',ncorn);
    elseif dim == 3
      C = textscan(fid,'%d %f64 %f64 %f64 %d %d %d',ncorn);
    end
      
    n = numel(C);
    ni.corn = [C{1},C{n-2:n}];
    xyz.cornx = C{2};
    xyz.corny = C{3};
    if dim == 3, xyz.cornz = C{4}; end
    
    %  Get side
    while isempty(strfind(str,'SIDES'))
      str = fgetl(fid);
    end
    nside = sscanf(str,'%g');
    
    %  Side index, a packed flag word (1=warped; 2=periodic; 4=contact;
    %  16=joint; 32=edge; 48=face; 64=exterior), the boundary identifier,
    %  the coefficient indices associated with this side (varying numbers,
    %  depending on dimension and order)
    if dim == 2
      C = textscan(fid,'%d %d %d %d',nside);
    elseif dim == 3
      C = textscan(fid,'%d %d %d',nside);
    end
    ni.side = [C{:}];
    
    %  Get legs (only for 3d)
    if dim == 3
      while isempty(strfind(str,'LEGS'))
        str = fgetl(fid);
      end
      nlegs = sscanf(str,'%g');
      C = textscan(fid,'%d %d %d',nlegs);
      ni.legs = [C{:}];
    end
    
    %  Get cell info (implicitly includes connectivity info
    while isempty(strfind(str,'CELLS'))
      str = fgetl(fid);
    end
    ncell = sscanf(str,'%g');
    
    %  Cell index, a packed flag word  (1=bounding; 2=warped), three corner
    %  /joint node indices, three side node indices, neighboring cell
    %  indices (correspond to listed sides...0 means no neighboring cell),
    %  region number and material number...note coefficient indices
    %  mentioned in documentation (they're supposed to be there) but they
    %  don't seem to appear in the data file
    %  Old way
%     C = textscan(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d',ncell);
    %  Better (we think), new way
    str = fgetl(fid);
    tmp = sscanf(str,'%d');
    nentries_per_line = numel(tmp);
    tmpstr = repmat('%d ',1,nentries_per_line); tmpstr = tmpstr(1:end-1);
    C = textscan(fid,tmpstr,ncell-1);
    cellinfo = [tmp';C{:}];
    
    %  Expand as need be in the future
    switch dim
      case 2 % Two dimensions
        cornpercell = 3;
        sidepercell = dof-cornpercell;
        switch order
          case 2 % Quadratic
            subcellpercell = 4;
            %  Cell info list comes with corner nodes then side nodes.
            %  This gives the subcells for that list of, in this case, 6
            %  nodes
            cellinds = [1,2,3];
            subcellinds = [1,6,5;6,2,4;6,4,5;4,3,5];
        end
      case 3
        %  Assume linear elements
        cornpercell = 4;
        sidepercell = 0;
        cellinds = [1,2,3,4];
        subcellinds = [1,2,3,4];
    end
    
    cellcorn = cellinfo(:,3:2+cornpercell);
    cellside = cellinfo(:,3+cornpercell:2+dof);
    
    %  Get cell info
    for rc = 1:nreg
      mexit = false;
      while isempty(strfind(str,'COEFFICIENTS')) && (~mexit)
        str = fgetl(fid);
        %  Special case.  There are no variables.  Leave.  Only good for
        %  linear meshes, I believe, as some side node stuff is done below,
        %  I think.
        if ~isempty(strfind(str,'ENDFILE'))
          mexit = true;
        end
      end
      if mexit
        nnode = sscanf(str,'%g');
      end
      
      %  "Coefficient" index, coefficients.  I assumed the coefficient was
      %  defined in the usual way with u(x_i) = u_i at the ith node but this
      %  seems to be wrong.  See documentation.
      mystr = '%d';
      for n = 1:nvar
        mystr = [mystr,' %64f'];
        if n < nnode, mystr = [mystr,' ']; end
      end
      C = textscan(fid,mystr,nnode);
      coefinfo(C{1},:) = [C{2:end}];
      str = 'Keep going';
    end
    if ~isempty(coefinfo)
      coefcorn = coefinfo(ni.corn(:,end),:);
      %  Adjust coefs on sides
      coefsides = zeros(nside,nvar);
    else
      coefcorn = []; coefsides = [];
    end
      
    %  Again, expand as need be

    switch order
      case 2
        %  Note some redundancy, many side values are calculated/stored
        %  more than once (any shared sides between elements, basically)
        ni.all = zeros(ncorn+nside,1);
        for c = coords
          xyz.(c) = zeros(ncorn+nside,1);
          xyz.(['side',c])(cellside(:,1)) = ...
            mean(xyz.(['corn',c])(cellcorn(:,[2,3])),2);
          xyz.(['side',c])(cellside(:,2)) = ...
            mean(xyz.(['corn',c])(cellcorn(:,[3,1])),2);
          xyz.(['side',c])(cellside(:,3)) = ...
            mean(xyz.(['corn',c])(cellcorn(:,[1,2])),2);
        end
        coefsides(cellside(:,1),:) = ...
          (coefcorn(cellcorn(:,2),:)+coefcorn(cellcorn(:,3),:))./2;
        coefsides(cellside(:,2),:) = ...
          (coefcorn(cellcorn(:,3),:)+coefcorn(cellcorn(:,1),:))./2;
        coefsides(cellside(:,3),:) = ...
          (coefcorn(cellcorn(:,1),:)+coefcorn(cellcorn(:,2),:))./2;
        coefsides = coefsides+coefinfo(ni.side(:,end),:);
        ni.all(ni.corn(:,end)) = ni.corn(:,2)>0;
        ni.all(ni.side(:,end)) = ni.side(:,2)>0;
        for c = coords
          xyz.(c)(ni.corn(:,end)) = xyz.(['corn',c]);
          xyz.(c)(ni.side(:,end)) = xyz.(['side',c]);
        end
    end
    if ~isempty(coefinfo)
      coefinfo(ni.side(:,end),:) = coefsides;
    end
    
%     %  For debugging
%     close(figure(1)); figure(1);
%     for rc = 1:size(cellcorn,1)
%       xc = xyz.cornx(cellcorn(rc,:));
%       xs = [xc(1)+xc(2),xc(2)+xc(3),xc(3)+xc(1)]./2;
%       yc = xyz.corny(cellcorn(rc,:));
%       ys = [yc(1)+yc(2),yc(2)+yc(3),yc(3)+yc(1)]./2;
%       plot(xc,yc,'b*-');
%       hold on
%       plot(xs,ys,'rs');
%       hold off
%       keyboard
%     end
    
    origcell = cellinfo(:,3:2+dof);
    for cc = 1:cornpercell
      origcell(:,cc) = ni.corn(cellinfo(:,2+cc),end);
    end
    for cc = cornpercell+[1:sidepercell]
      origcell(:,cc) = ni.side(cellinfo(:,2+cc),end);
    end
    
    origtetra = reshape(origcell(:,cellinds')',numel(cellinds),[])';
    tetra = reshape(origcell(:,subcellinds')',numel(subcellinds),[])';
    
    %  Find the indices associated with each tetrahedra.  Put them in the
    %  ni.side structure
    if dim == 3
      faceinds = [1 2 3;1 4 2;2 4 3;1 3 4];
      faceinds = faceinds([1 4 2 3],:);
      for cc = 1:ncell
        for fc = 1:4
          sideind = cellinfo(cc,6+fc);
          ni.side(sideind,4:6) = cellinfo(cc,2+faceinds(fc,:));
        end
      end
    end
    
%     %  For debugging purposes only
%     figure(1); close(figure(1));
%     for tc = 1:size(tetra,1)
%       plot(xyz.x(tetra(tc,:)),xyz.y(tetra(tc,:)));
%       hold on
%       pause
%     end
    
    for vc = 1:nvar
      rest.(var_name{vc}) = coefinfo(:,vc);
    end
    rest.reg_name = regs{rc};
  end
  
%   if nreg > 1
%     xyz = xyzs;
%     rest = rests;
%     DT = DTs;
%     tetra = tetras;
%     tri = tris;
%     t0 = t0s;
%     origtetra = origtetras;
%     node_info = node_infos;
%   else
%     xyz = xyzss{1};
%     rest = rests{1};
%     DT = DTs{1};
%     tetra = tetras{1};
%     tri = tris{1};
%     t0 = t0s{1};
%     origtetra = origtetras{1};
%     node_info = node_infos{1};
%   end
    
  fclose(fid);

end
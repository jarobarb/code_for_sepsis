function [myhandle,myfield] = tetra_exterior_plot(currfield,xyz,tri,...
	node_info,varargin)
  
  %  Typical usage:
	%  [myhandle,myfield] = tetra_exterior_plot(currfield,xyz,tri,...
	%    node_info)
	%  To plot material 1 using data from read_in_fpde_dat:
% 	 [myhandle,myfield] = tetra_exterior_plot(rest{1}.U,xyzs{1},tri{1},...
% 	   node_info{1});
	%  currfield-vector containing field that one wishes to look at
	%  xyz-matrix with rows of xyz coords for nodes corresponding to
	%    currfield
	%  node_info-raw field data (see read_in_fpde_dat.m)
	%  myhandle-graphics handle so that one can adjust the patch's properties
	%    easily
	%  myfield-only the values of the current field at locations within x,y,
	%    and z bounds (see xbnds, ybnds, zbnds below)
	%
	%  Useful things to know
	%    colorbar puts a colorbar up
	%    caxis allows you to change the lower and upper limits of that colorbar
	%    myhandle is the graphics handle
	
  %  If tri is actually a bunch of tetrahedra...convert it into triangles
  if size(tri,2) == 4
    tetra = tri;
    tetratotri = [1,2,3;1,4,3;1,2,4;2,3,4];
    tri = reshape(tetra(:,tetratotri')',3,[])';
  end
  
  %  Allow one to limit the extent of the mesh in x/y/z directions
  xbnds = [];
  ybnds = [];
  zbnds = [];
  
  if numel(varargin) > 0
    xbnds = varargin{1};
  end
  if numel(varargin) > 1
    ybnds = varargin{2};
  end
  if numel(varargin) > 2
    zbnds = varargin{3};
  end
  
  if isfield(xyz,'z'), dim = 3; else dim = 2; end
  
  tmp = node_info(:,2);
  extinds = all(tmp(tri),2);
  if ~isempty(xbnds)
    extinds = extinds & all(xyz.x(tri) >= xbnds(1),2) & ...
      all(xyz.x(tri) <= xbnds(2),2);
  end
  if ~isempty(ybnds)
    extinds = extinds & all(xyz.y(tri) >= ybnds(1),2) & ...
      all(xyz.y(tri) <= ybnds(2),2);
  end
  if ~isempty(zbnds)
    extinds = extinds & all(xyz.z(tri) >= zbnds(1),2) & ...
      all(xyz.z(tri) <= zbnds(2),2);
  end
  trinew = tri(extinds,:)';
  myfield = currfield(trinew);
  if dim == 3
    myhandle = patch(xyz.x(trinew),xyz.y(trinew),xyz.z(trinew),currfield(trinew));
  else
    myhandle{1} = patch(xyz.x(tri'),xyz.y(tri'),currfield(tri'));
    extinds = find(tmp > 0);
    ec_ind = extinds(1);
    myinds = extinds(1)*ones(size(extinds));
    for ec = 2:numel(extinds)
      tmp2 = find(any(tri == ec_ind,2));
      tmp3 = tri(tmp2,:);
      neighbor_nodes = unique(tmp3(:));
      neighbor_edge_nodes = intersect(setdiff(neighbor_nodes,ec_ind),...
        extinds);
      if ~ismember(neighbor_edge_nodes(1),myinds)
        myinds(ec) = neighbor_edge_nodes(1);
      else
        myinds(ec) = neighbor_edge_nodes(2);
      end
      ec_ind = myinds(ec);
    end
    myinds = [myinds;myinds(1)];
    hold on
    myhandle{2} = plot(xyz.x(myinds),xyz.y(myinds),'k');
  end
  axis tight
  axis equal
  set(gcf,'Renderer','zbuffer');

end
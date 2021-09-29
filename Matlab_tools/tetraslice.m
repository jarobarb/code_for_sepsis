function [hs,xyzout,valorvec] = tetraslice(xyzin,valorvec,tetra,normal,pt,varargin)

  %  [hs,xyzout,valorvec] = tetraslice(xyzin,valorvec,tetra,normal,pt,...
  %    varargin)
  %
	%  Typical usage (Stress_2d.dat):
	%  [hs,xyzout,valorvec] = tetraslice(xyzs{1},rest{1}.U,tetra{1},...
	%    [0,0,1],[0,0,0]);
	%  [hs,xyzout,valorvec] = tetraslice(xyzs{1},rest{1}.U,tetra{1},...
	%    [-1,1,0],[0,0,0]);
	%  Notice, using both commands one right after the other will make it so
	%    that both slices stay on the plot.
	%
  %  Given some information about a tetrahedral mesh and values of
  %  variables defined on that mesh, tetraslice allows one to draw contour
  %  plots of the variables on planes that slice through the tetrahedral
  %  mesh.  One can also project on to cylinders using parameter pairs
  %  'r',radiusval as well.
  %
  %  OUTPUTS
  %  hs-cell of graphics handles for use in adjusting the plots (e.g.
  %    colormap, edgecolor, etc.).  For contours there is one handle for
  %    triangles, another for rectangles, another for pentagons, etc
  %    (depending on how many of each shape are obtained when the
  %    tetrahedrons are intersected with the plane
  %  xyzout-matrix [x,y,z] of all the intersection points of the tetrahedra
  %    with the planes
  %  valorvec-column vector or matrix of the corresponding values of the
  %    variables at the intersection points of the tetrahedra with the
  %    planes
  %  
  %  INPUTS
  %  xyzin-matrix of x, y, z components of a set of vertices for a
  %    tetrahedral mesh [x,y,z] or a structure with fields .x, .y, .z that
  %    give x, y, z components of the set of vertices.
  %  valorvec-either 1. a column vector containing the values of a single
  %    variable defined on the mesh for use in making a contour plot, 2. a
  %    matrix of [u,v,w] where u, v, and w are the x, y, and z components
  %    of a particular vector (e.g. velocity, current, heat flux) that you
  %    want to plot, 3.  a structure with some fields that correspond to
  %    the variable(s) you want to use the program to plot.  If 3. then one
  %    must also hand in a parameter value pair
  %    'field_names',{'field1'} for contour plots or the parameter value
  %    pair 'field_names',{'field1','field2','field3'} for vector plots
  %  tetra-nx6 matrix where each row has the vertex indices for each
  %    tetrahedron...that is x(tetra(i,:)) gives the x values for the
  %    vertices of the 6 points on the ith tetrahedron
  %  normal-a normal vector to the plane of interest
  %    for cylinders, vector parallel to cylinder axis
  %  pt-a point in the plane of interest
  %    for cylinders, point on cylinder axis
  %  varargin-allows one to feed in additional options like the use of the
  %    field_names option mentioned above.  See m-file for more options
  
  %  Default options
  %  One doesn't even have to make a plot...this can be used just to return
  %  estimates for the field or fields at various points along the plane of
  %  interest.
  do_plot = true;
  %  Whether or not to put on the colorbar
  do_colorbar = true;
  %  This allows one to hand in a structure for valorvec and then use the
  %  varargin pair 'field_names',{'field1'} for contour plots of field1
  %  or the varargin pair 'field_names',{'field1','field2','field3'} for
  %  vector plots using the vectors [field1,field2,field3].
  field_names = {};
  %  Replaces default values if appropriate parameter pair is defined via
  %  varargin
  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  %  Allows one to hand in a structure instead of a matrix for the xyz info
  if isstruct(xyzin)
    xyzin = [xyzin.x(:),xyzin.y(:),xyzin.z(:)];
  end
  %  Converts valorvec from a structure to a column vector or matrix as
  %  appropriate...must define field_names via parameter pairs in order for
  %  this to work (see above)
  valorvecnew = [];
  if ~isempty(field_names)
    for fnc = 1:numel(field_names)
      valorvecnew = [valorvecnew,valorvec.(field_names{fnc})];
    end
  end
  if ~isempty(valorvecnew)
    valorvec = valorvecnew;
  end
  
  %  Allows user to input a radius "r"...if this is input then we assume
  %  that we want to project the solution onto a cylinder rather than a
  %  plane.  This makes "normal" the vector tangent to the axis of the
  %  cylinder and "pt" a point the cylinder axis passes through  
  dist_to_plane = (normal(1)*(xyzin(:,1)-pt(1))+...
    normal(2)*(xyzin(:,2)-pt(2))+...
    normal(3)*(xyzin(:,3)-pt(3)))./norm(normal);
  if exist('r')    
    hypoteneuse2 = (xyzin(:,1)-pt(1)).^2+(xyzin(:,2)-pt(2)).^2+...
      (xyzin(:,3)-pt(3)).^2;
    dist_to_line = sqrt(hypoteneuse2-dist_to_plane.^2);
    %  dist_to_plane now effectively becomes distance to the cylindrical
    %  surface
    dist_to_plane = dist_to_line-r;
  end

  %  Second set of tetrahedron indices including only tetrahedron that
  %    straddle the plane/cylinder.  We call these straddling tetrahedra.
  tetra2 = tetra(any((dist_to_plane(tetra) >= 0)')' & ...
    any((dist_to_plane(tetra) <= 0)')',:);
  
  %  Pick out just the distances of the vertices of the straddling
  %    tetrahedra to the plane/cylinder.
  dist_to_plane_tetra2 = dist_to_plane(tetra2);
  
  %  Each straddling tetrahedron has 6 edges.  Some of those will intersect
  %  the plane while others may not.  We initialize matrices that can be
  %  used to find those intersecting edges and the x, y, z and variable
  %  values at those edges.
  edge_intersection(size(dist_to_plane_tetra2,1),6) = false;
  intersection_x(size(dist_to_plane_tetra2,1),6) = 0;
  intersection_y(size(dist_to_plane_tetra2,1),6) = 0;
  intersection_z(size(dist_to_plane_tetra2,1),6) = 0;
  for fnc = 1:size(valorvec,2)
    intersection_val{fnc}(size(dist_to_plane_tetra2,1),6) = 0;
  end
  c3 = 0;
  for c1 = 1:3
    for c2 = c1+1:4
      c3 = c3+1;
      edge_intersection(:,c3) = ...
        (dist_to_plane_tetra2(:,c1).*dist_to_plane_tetra2(:,c2)<=0) &...
        ~((dist_to_plane_tetra2(:,c1) == 0) & ...
        (dist_to_plane_tetra2(:,c2) == 0));
      intersection_x(:,c3) = (dist_to_plane_tetra2(:,c1).*...
        xyzin(tetra2(:,c2),1)-dist_to_plane_tetra2(:,c2).*...
        xyzin(tetra2(:,c1),1))./(dist_to_plane_tetra2(:,c1)-...
        dist_to_plane_tetra2(:,c2));
      intersection_y(:,c3) = (dist_to_plane_tetra2(:,c1).*...
        xyzin(tetra2(:,c2),2)-dist_to_plane_tetra2(:,c2).*...
        xyzin(tetra2(:,c1),2))./(dist_to_plane_tetra2(:,c1)-...
        dist_to_plane_tetra2(:,c2));
      intersection_z(:,c3) = (dist_to_plane_tetra2(:,c1).*...
        xyzin(tetra2(:,c2),3)-dist_to_plane_tetra2(:,c2).*...
        xyzin(tetra2(:,c1),3))./(dist_to_plane_tetra2(:,c1)-...
        dist_to_plane_tetra2(:,c2));
      for fnc = 1:size(valorvec,2)
        intersection_val{fnc}(:,c3) = (dist_to_plane_tetra2(:,c1).*...
          valorvec(tetra2(:,c2),fnc)-dist_to_plane_tetra2(:,c2).*...
          valorvec(tetra2(:,c1),fnc))./(dist_to_plane_tetra2(:,c1)-...
          dist_to_plane_tetra2(:,c2));
      end
    end
  end
  
  %  Store the intersection locations and values at those locations
  %  Note there is no organization whatsoever to these...all tetrahedral
  %  info is lost and these become random points with random values located
  %  at those points.
  %  Note valorvec not used after this point...hence the redefinition is
  %  currently ok.
  x = intersection_x(edge_intersection);
  y = intersection_y(edge_intersection);
  z = intersection_z(edge_intersection);
  xyzout = [x(:),y(:),z(:)];
  valorvec = [];
  for fnc = 1:numel(intersection_val)
    tmp = intersection_val{fnc}(edge_intersection);
    valorvec = [valorvec,tmp(:)];
  end
  
  %  Note hs is defined as a matrix rather than a cell...this is because I
  %  believe I have other programs that utilize this difference.
  if ~do_plot, hs = []; return; end
  
  %  Quadrilaterals and pentagons can get twisted so we temporarily project
  %  the vertices onto a 2d space (u,v space here), use convhull to untwist
  %  them, then store them in a matrix/list of untwisted polygons all with
  %  the same number of vertices.
  uvvec = null(normal);
  u = uvvec(:,1);
  v = uvvec(:,2);
  
  if size(valorvec,2) > 1
%     hs = quiver3(x,y,z,valorvec(:,1),valorvec(:,2),valorvec(:,3));
    us = u(1)*x+u(2)*y+u(3)*z;
    vs = v(1)*x+v(2)*y+v(3)*z;
    fu = u(1)*valorvec(:,1)+u(2)*valorvec(:,2)+u(3)*valorvec(:,3);
    fv = v(1)*valorvec(:,1)+v(2)*valorvec(:,2)+v(3)*valorvec(:,3);
    hs = jared_quiver(us,vs,fu,fv,'rescale',10);
    axis tight
    axis equal
  end
  
  %  Each straddling tetrahedron can intersect the plane/cylinder in
  %  various ways forming a triangle, quadrilateral, or maybe a pentagon
  %  (though I think it is hard to find pentagons and impossible to find
  %  heaxagons).  To vectorize the patch command, we plot each shape
  %  separately.  Below gives the number of vertices for each straddling
  %  tetrahedron and tells us, for each straddling tetrahedron, if the
  %  tetrahedron intersects in a shape with 3 vertices (triangle), 4
  %  vertices (quadrilateral) or more
  vertices_per_poly = unique(sum(edge_intersection,2));
  
  if numel(intersection_val) == 1
    intersection_val = intersection_val{1};
    for vc = vertices_per_poly'
      
      inds = bsxfun(@and,sum(edge_intersection,2) == vc,edge_intersection)';
      
      edge_intersection = edge_intersection';
      intersection_x = intersection_x';
      intersection_y = intersection_y';
      intersection_z = intersection_z';
      intersection_val = intersection_val';
      
      tmpx = reshape(intersection_x(inds),vc,[]);
      tmpy = reshape(intersection_y(inds),vc,[]);
      tmpz = reshape(intersection_z(inds),vc,[]);
      tmpval = reshape(intersection_val(inds),vc,[]);
      
      if vc > 3
        if ~exist('r')
          tmpu = (tmpx-pt(1))*u(1)+(tmpy-pt(2))*u(2)+(tmpz-pt(3))*u(3);
          tmpv = (tmpx-pt(1))*v(1)+(tmpy-pt(2))*v(2)+(tmpz-pt(3))*v(3);
        else
          %  If a cylinder, u becomes along the axial direction while v is
          %  theta
          tmpu = (tmpx-pt(1))*u(1)+(tmpy-pt(2))*u(2)+(tmpz-pt(3))*u(3);
          tmpv = (tmpx-pt(1))*v(1)+(tmpy-pt(2))*v(2)+(tmpz-pt(3))*v(3);
          tmpz2 = (tmpx-pt(1))*normal(1)+(tmpy-pt(2))*normal(2)+...
            (tmpz-pt(3))*normal(3);
          tmptheta = atan2(tmpu,tmpv);
          tmpu = tmpz2;
          tmpv = tmptheta;
        end
        tmpa(size(tmpu,2)) = 0;
        tmplist = true(1,size(tmpu,2));
        
        for qc = 1:size(tmpu,2)
          try
            [K,tmpa(qc)] = convhull(tmpu(:,qc),tmpv(:,qc));
            if numel(K) < vc
              K(end+1:end+vc) = K(1);
            end
          catch
            tmpa(qc) = 0;
          end
          
          if tmpa(qc) > 1e-15
            tmpu(:,qc) = tmpu(K(1:vc),qc);
            tmpv(:,qc) = tmpv(K(1:vc),qc);
            tmpx(:,qc) = tmpx(K(1:vc),qc);
            tmpy(:,qc) = tmpy(K(1:vc),qc);
            tmpz(:,qc) = tmpz(K(1:vc),qc);
            tmpval(:,qc) = tmpval(K(1:vc),qc);
          else
            tmplist(qc) = false;
          end
        end
      end
      current_h_ind = find(vc == vertices_per_poly);
      
      hs{current_h_ind} = patch(tmpx,tmpy,tmpz,tmpval);
      set(hs{current_h_ind},'EdgeColor','none');%'interp');
      
      if vc < vertices_per_poly(end)
        edge_intersection = edge_intersection';
        intersection_x = intersection_x';
        intersection_y = intersection_y';
        intersection_z = intersection_z';
        intersection_val = intersection_val';
      end
      
    end
    
    axis equal
    axis tight
      
    if do_colorbar
      colorbar
    end
      
    view(-normal)
    
  end    
  
end
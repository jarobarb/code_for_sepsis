function [junk1,vertices,y,junk2,junk3] = ...
  par_simplex2(vertices,fhereiam,varargin)

  %  Junk variables for compatibility with par_simplex
  junk1 = []; junk2 = []; junk3 = [];

  %  Parallel algorithm taken from Lee and Wiswall 2007:  "A parallel
  %  implementation of the simplex function minimization routine"

  %  Default values
  tol = 1e-3;
  parcomp = true;
  maxits = 1e6;
  bump = 0.01;

  for vac = 1:2:numel(varargin)
    eval([varargin{vac},' = varargin{vac+1};']);
  end
  
  %  If we are inside another parallelized process, don't try to
  %  parallelize this process, Matlab will fail if we do.
  if isinparallel, parcomp = false; end
  parcomp = false;

  %  Make a simplex if only one vertex is handed in
  if min(size(vertices)) == 1
    vertices = vertices(:)';
    vertices = [vertices;
      vertices.*(ones(numel(vertices),1)+bump*eye(numel(vertices)))];
  end
  %  Number of free parameters
  ndim = min(size(vertices));

  %  Evaluate initial simplex
  if parcomp
    p = gcp();
    nproctouse = min(p.NumWorkers,size(vertices,2));
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
    nproctouse = 1;
    for n = 1:ndim+1
      y(n) = fhereiam(vertices(n,:)');
      fprintf('%g %d/%d done.\n',y(n));
    end
  end
  
  vertex_dist_bound = max(pdist(vertices));

  while vertex_dist_bound > tol
    %  Some nice monitors
    fprintf('Rel max vert dist bnd:  %18.16g.  Min y values:  ',...
      vertex_dist_bound/tol);
    fprintf('%18.16g %18.16g',mink(y,2)); fprintf('\n');
    
    %  Resort so the min val is at the front of y/top of the vertices and
    %  the max val is at the back of y/bottom of the vertices
    [y,inds] = sort(y);
    vertices = vertices(inds,:);
    
    %  Store vertices
    vertices_new = vertices;
    y_new = y;
    
    %  Parfors may not be the best route.  See above for alternate route
    %  which may take less time.
    centroid = mean(vertices(1:end-nproctouse,:),1);
    parfor pc = 1:nproctouse
      %  Use of y(end-nproctouse) may be different than the Lee+Wiswall
      %  algorithm (they use y(end-pc)) but what is below seemed to produce
      %  better results quicker for our particular test case (see
      %  test_par_simplex2)***changed back 8/18/2020
      %  Calculate the centroid
%       inds = [1:(ndim+1)]~=(ndim+1-pc+1);
%       centroid(pc,:) = mean(vertices(inds,:),1);
%     
%       [tmp1{pc},tmp2{pc}] = ...
%         flipper(fhereiam,centroid(pc,:),vertices(end-pc+1,:),...
%         y(end-pc+1),y(1),y(end-nproctouse));
      
      [tmp1{pc},tmp2{pc}] = ...
        flipper(fhereiam,centroid,vertices(end-pc+1,:),...
        y(end-pc+1),y(1),y(end-pc+1));
    end
    for pc = 1:nproctouse
      vertices_new(end-pc+1,:) = tmp1{pc};
      y_new(end-pc+1) = tmp2{pc};
    end
    
    %  Nothing changed...shrink towards best point
    if isequal(vertices_new,vertices)
      vertices(2:end,:) = (vertices(2:end,:)+vertices(1,:))./2;
      parfor vc = 1:size(vertices(2:end,:),1)
        y(vc) = fhereiam(vertices(vc,:)');
      end
    else
      vertices = vertices_new;
      y = y_new;
    end
        
    vertex_dist_bound = max(pdist(vertices));
  end

end

function [vertex_new,y_new] = ...
  flipper(f,centroid,vertex_old,y_old,y_best,y_near_worst)

  reflected_vertex = centroid+1*(centroid-vertex_old);
  y_reflected = f(reflected_vertex');
  
  if y_reflected < y_best
    reflected_expanded_vertex = reflected_vertex+1*(reflected_vertex-...
      centroid);
    y_reflected_expanded = f(reflected_expanded_vertex');
    if y_reflected_expanded < y_reflected
      vertex_new = reflected_expanded_vertex;
      y_new = y_reflected_expanded;
    else
      vertex_new = reflected_vertex;
      y_new = y_reflected;
    end
  elseif y_reflected < y_near_worst
    vertex_new = reflected_vertex;
    y_new = y_reflected;
  else
    if y_reflected < y_old
      y_best_so_far = y_reflected;
      vertex_best_so_far = reflected_vertex;
      contracted_vertex = (centroid+reflected_vertex)/2;
    else
      y_best_so_far = y_old;
      vertex_best_so_far = vertex_old;
      contracted_vertex = (centroid+vertex_old)/2;
    end
    y_contracted = f(contracted_vertex');
    if y_contracted < y_best_so_far
      y_new = y_contracted;
      vertex_new = contracted_vertex;
    else
      y_new = y_best_so_far;
      vertex_new = vertex_best_so_far;
    end
  end
 
end

function answer = isinparallel()

  try
    answer = ~isempty(getCurrentTask());
  catch err
    if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
      rethrow(err);
    end
    answer = false;
  end
  
end
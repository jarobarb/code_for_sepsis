myf = @(x) sum(x.^2);
[vertices,y] = par_simplex2((1:10)/10,myf,'parcomp',true,'tol',1e-6,...
  'bump',0.01);
[vertices,y] = par_simplex2(vertices(1,:),myf,'parcomp',true,'tol',1e-6,...
  'bump',0.01);

% [a1,a2,a3,a4,a5] = par_simplex((1:10)/10,myf);
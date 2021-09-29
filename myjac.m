function J = myjac(t,y,rp,op,fp,strp)

  %  1-B-Bacteria, 2-M-Macrophages/pro-inflammatory, 
  %  3-A-anti-inflammatory, 4-epsilon-damage
%   y = 0.1*rand(size(y)); y(1) = 2*fp.clot_capacity*rand;
%   y(5) = 2*fp.clot_capacity*rand;
  B = y(1);
  M = y(2);
  A = y(3);
  epsilon = y(4);
  Bc = y(5);
  
  eq2a = rp.kM*M+rp.kB*op.c1*B/fp.Bmf+rp.kepsilon*epsilon;
  eq3a = y(2)+rp.kAMEps*y(4);
  est_flux = (fp.kD1.*y(5)-fp.kD2.*y(1));
  flux_ind = (~fp.limit_flux) || ((y(5) < fp.clot_capacity) ||...
    (est_flux > 0));
  clot_to_body_flux = est_flux.*flux_ind;
  carrying_capacity = fp.effective_carrying_capacity*fp.B_infy+...
    (~fp.effective_carrying_capacity)*fp.B_infy*fp.Bmf;
  
  J = fp.tscale*[op.k1*(1-2*B/carrying_capacity)-...
    (rp.kBl*rp.sl)*fp.Bmf^2*rp.mul/(B*rp.klB+fp.Bmf*rp.mul)^2-...
    fp.k5*M/(A*fp.kA+1)-flux_ind*fp.kD2,...
    -fp.k5*B/(A*fp.kA+1),...
    fp.k5*B*M*fp.kA/(A*fp.kA+1)^2,...
    0,...
    fp.kD1.*flux_ind;
    
    rp.nu1*rp.nu2*rp.kB*op.c1/((A*fp.kA+1)*(rp.nu2+eq2a)^2*fp.Bmf),...
    rp.nu1*rp.nu2*rp.kM/((A*fp.kA+1)*(rp.nu2+eq2a)^2)-rp.muM,...
    -rp.nu1*eq2a*fp.kA/((A*fp.kA+1)^2*(rp.nu2+eq2a)),...
    rp.nu1*rp.nu2*rp.kepsilon/((A*fp.kA+1)*(rp.nu2+eq2a)^2),...
    0;
    
    0,...
    rp.a1*fp.nu4/((A*fp.kA+1)*(fp.nu4+eq3a)^2),...
    -rp.a1*fp.kA*eq3a/((A*fp.kA+1)^2*(fp.nu4+eq3a))-rp.muA,...
    rp.a1*fp.nu4*rp.kAMEps/((A*fp.kA+1)*(fp.nu4+eq3a)^2),...
    0;
    
    0,...
    (fp.f*M>op.T)*fp.f/(A*fp.kA+1),...
    -max(fp.f*M-op.T,0)*fp.kA/(A*fp.kA+1)^2,...
    -1/op.tau,...
    0;
    
    fp.kD2.*flux_ind,...
    0,...
    0,...
    0,...
    -fp.kD1.*flux_ind+fp.kD3*(1-2*fp.use_log_Bc*Bc/fp.clot_capacity)];
  
%   J2 = J;
%   
%   baseline = dydt(t,y,rp,op,fp,strp);
%   for vc = 1:numel(y)
%     ytmp = y;  ytmp(vc) = ytmp(vc).*(1+sqrt(eps));
%     new = dydt(t,ytmp,rp,op,fp,strp);
%     J2(:,vc) = (new-baseline)./(ytmp(vc)-y(vc));
%   end
%   
%   [fp.limit_flux, flux_ind, fp.use_log_Bc]
%   [J;J2;(J-J2)./(abs(J)+eps(J))]

end
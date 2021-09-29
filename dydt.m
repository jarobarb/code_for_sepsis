function [ftys,Js] = dydt(t,ys,rps,ops,fps,strps)

%     y(1) = fp.Bscale*y(1); y(2) = fp.Mscale*y(2);
%     y(3) = fp.Ascale*y(3); y(4) = fp.escale*y(4);
%     y = max(y,0);
    %  The right hand side (i.e. f(t,\vec{y})) of the system of diff eqs.
    
    ftys = zeros(size(ys));
    Js = zeros(numel(ys));
%     Js = sparse(numel(ys),numel(ys));
    
    for cc = 1:numel(rps)
      rp = rps(cc); op = ops(cc); fp = fps(cc); strp = strps(cc);
      inds = 4*(cc-1)+[1:5];
      y = ys(inds);
      
      eq2a = rp.kM*y(2)+rp.kB*op.c1*y(1)/fp.Bmf+rp.kepsilon*y(4);
      eq3a = y(2)+rp.kAMEps*y(4);
      est_flux = (fp.kD1.*y(5)-fp.kD2.*y(1));
      flux_ind = (~fp.limit_flux) || ((y(5) < fp.clot_capacity) ||...
        (est_flux > 0));
      clot_to_body_flux = est_flux.*flux_ind;
      carrying_capacity = fp.effective_carrying_capacity*fp.B_infy+...
        (~fp.effective_carrying_capacity)*fp.B_infy*fp.Bmf;
      
      ftys(inds) = fp.tscale.*...
        [op.k1*y(1)*(1-y(1)/(carrying_capacity))-...
        (rp.kBl*rp.sl)*y(1)/(rp.mul+rp.klB*(y(1)/fp.Bmf))-...
        fp.k5*y(1)*y(2)/(1+fp.kA*y(3))+...
        clot_to_body_flux;
        
        fp.sM+rp.nu1*eq2a/(1+fp.kA*y(3))/(rp.nu2+eq2a)-rp.muM*y(2);
        
        rp.sA+rp.a1*eq3a/(1+fp.kA*y(3))/(fp.nu4+eq3a)-rp.muA*y(3);
        
        (op.epsilon0-y(4))/op.tau+max(fp.f*y(2)-op.T,0)/(1+fp.kA*y(3));...
        -clot_to_body_flux+fp.kD3*y(5)*...
        (1-fp.use_log_Bc*y(5)/fp.clot_capacity)];
      
      if nargout > 1
        Js(inds,inds) = myjac(t,ys(inds),rp,op,fp,strp);
      end
    end

end
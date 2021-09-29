function [fty,J] = dydt_for_eq3(t,y,rp,op,fp,strp)

%     y(1) = fp.Bscale*y(1); y(2) = fp.Mscale*y(2);
%     y(3) = fp.Ascale*y(3); y(4) = fp.escale*y(4);
%     y = max(y,0);
    %  The right hand side (i.e. f(t,\vec{y})) of the system of diff eqs.
    
    B = y(1); M = y(2); Bc = y(3);
    
    k1 = op.k1;
    Bs = 1/fp.Bmf;
    Binfty = fp.B_infy;
    k2 = rp.kBl;
    sl = rp.sl;
    k3 = rp.klB;
    mul = rp.mul;
    k5 = fp.k5;
    nu1 = rp.nu1;
    kB = rp.kB;
    c1 = op.c1;
    kM = rp.kM;
    keps = rp.kepsilon;
    nu2 = rp.nu2;
    muM = rp.muM;
    muA = rp.muA;
    sA = rp.sA;
    a1 = rp.a1;
    k4 = rp.kAMEps;
    tau = op.tau;
    f = fp.f;
    T = op.T;
    kA = fp.kA;
    kD1 = fp.kD1;
    kD2 = fp.kD2;
    kD3 = fp.kD3;
    Cc = fp.clot_capacity;
    %  gamma = 0:  Isler background anti-inflammatory levels = 0
    %  gamma = 1:  Isler background anti-inflammatory levels = in vivo
    %    background levels = model estimate of A_back = sA/muA
    %  gamma = 0.5:  Isler background anti-inflammatory levels = 1/2 in
    %    background levels = (sA/muA)/2
    %  We think gamma = 0.5 allows us to most comfortably say that
    %    interactions are "inhibited" 75% as 0.5 corresponds to an average
    %    of relative and absolute inhibition levels in our system, which we
    %    think agrees with Reynolds et al
    gamma = op.gamma;
    frac = (1-0.75)/(1+gamma*kA*sA/muA);
    
    est_flux = (kD1.*Bc-kD2.*B);
    flux_ind = (~fp.limit_flux) || ((Bc < fp.clot_capacity) ||...
      (est_flux > 0));
    clot_to_body_flux = est_flux.*flux_ind;
    use_log_Bc = fp.use_log_Bc;
    carrying_capacity = fp.effective_carrying_capacity*fp.B_infy+...
        (~fp.effective_carrying_capacity)*fp.B_infy*fp.Bmf;
        
    fty = ...
      [k1*B*(1-B/carrying_capacity)-k2*sl*B/(B*Bs*k3+mul)-k5*B*M*frac+...
      clot_to_body_flux;
      nu1*(kB*c1*Bs*B+kM*M+keps*(M*f-T)*frac*tau)*frac/(kB*c1*Bs*B+kM*M+...
      keps*(M*f-T)*frac*tau+nu2)-muM*M;
      -flux_ind*(-B*kD2+Bc*kD1)+kD3*Bc*(1-use_log_Bc*Bc/Cc)];
    
    J = [k1*(1-B*Bs/Binfty)-k1*B*Bs/Binfty-k2*sl/(B*Bs*k3+mul)+...
      k2*sl*B/(B*Bs*k3+mul)^2*Bs*k3-k5*M*frac-flux_ind*kD2,...
      -k5*B*frac,...
      flux_ind*kD1;
      
      nu1*kB*c1*Bs*frac*nu2/(kB*c1*Bs*B+kM*M+keps*(M*f-T)*frac*tau+...
      nu2)^2,...
      (-(kB*c1*Bs*B+kM*M+keps*(M*f-T)*frac*tau+nu2)^2*muM+...
      nu1*nu2*frac*(f*frac*keps*tau+kM))/(kB*c1*Bs*B+kM*M+keps*(M*f-T)*...
      frac*tau+nu2)^2,...
      0;
      
      flux_ind*kD2,...
      0,...
      ((-kD1*flux_ind+kD3)*Cc-2*kD3*Bc*use_log_Bc)/Cc];

end
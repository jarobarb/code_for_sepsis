function [fp,k,klabels,ysep] = solver_for_kA(rp,op,fp,strp,k,klabels,...
  Bmax,varargin)

  tol = 1e-15;
  blendfact = 1;
  percentage = op.inhi_perc;

  new_fp = @(kA) setfield(fp,'kA',kA);
  
%   myf = @(y) dydt_for_eq_kA(y,rp,op,new_fp(y(5)),strp);
%   
%   y0 = [Bmax,1,0.25,1,6.667];
%   
%   ysol = fsolve(@(y) myf(exp(y)),log(y0),...
%     optimoptions('fsolve','MaxFunctionEvaluations',1e5,...
%     'MaxIterations',1e5,...
%     'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,...
%     'Display','iter','algorithm','levenberg-marquardt'));
%   
%   y0 = [Bmax,1,1,1,6.667];
%   
%   myf = @(y) eq_sol_eqn(y,rp,op,fp,strp,0);
%   
%   ysol = fminsearch(@(y) sum((myf(y)).^2'+...
%     (y < [100,0.1,0.35,0.1,0.1])+...
%     (y > [10000,10,10,10,20])),...
%     y0,optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5));
%   
  %  Given a kA, this finds the corresponding septic solution (hopefully)
%   ysol = [Bmax,1,1,1,6.667];
%   
%   myf = @(y) eq_sol_eqn(y,rp,op,fp,strp,1,ysol(end));
%       
%   ysol = fminsearch(@(y) sum((myf(y)).^2),...
%     ysol,optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5,...
%     'TolX',1e-14,'TolFun',1e-14));
%   
%   kAs = linspace(1,20,201);
%   
%   for kAc = 1:numel(kAs)
%     ysol(end) = kAs(kAc);
%     myf = @(y) eq_sol_eqn(y,rp,op,fp,strp,1,ysol(end));
%       
%     ysol = fminsearch(@(y) sum((myf(y)).^2),...
%       ysol,optimset('Display','iter','MaxFunEvals',1e5,'MaxIter',1e5,...
%       'TolX',1e-14,'TolFun',1e-14));
%     
%     [a(kAc,:),b(kAc),c(kAc)] = eq_sol_eqn(ysol,rp,op,fp,strp,1,ysol(end));
%     d(kAc,:) = ysol;
%     
%     fprintf('%g/%g\n',kAc,numel(kAs));
%     
%   end
%   
%   n = 100;
%   for nc = 0:n
%     myf = @(y) eq_sol_eqn(y,rp,op,fp,strp,1-nc/n);
%     ysol = fsolve(@(y) myf(exp(y)),log(ysol),...
%       optimoptions('fsolve','MaxFunctionEvaluations',1e5,...
%       'MaxIterations',1e5,...
%       'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,...
%       'Display','iter','algorithm','levenberg-marquardt'));
%     ysol = exp(ysol);
%   end
%   
%   fminsearch(@(y) sum(myf(exp(y)).^2),log(ysol));
%   
%   myf = @(y) dydt_for_eq(0,y,rp,op,fp,strp);
%   
%   y0 = [Bmax,100,100,100];
%   
%   ysol = fsolve(@(y) myf(exp(y)),log(y0),...
%     optimoptions('fsolve','MaxFunctionEvaluations',1e6,...
%     'MaxIterations',1e6,...
%     'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,...
%     'Display','iter','algorithm','levenberg-marquardt'))
  
%   keyboard
%   
%   tstart = tic;
%   y0 = [Bmax,0,0.25,0,6.667];
%   tolvec=1e-12.*ones(size(y0)); tolvec(end) = 1e-16;
%   opts = odeset('InitialStep',0.1,'NonNegative',ones(size(y0)),...
%       'MaxStep',10000,'RelTol',0.01*min(tolvec),...
%       'AbsTol',0.01*min(tolvec),...
%       'Events',@(t,y,rp,op,fp,strp) myevent_augmented(t,y,rp,op,fp,strp,...
%       0.01*tolvec,tstart));
%     odes_sol = ode15s(@dydt_augmented,[0,1/eps],y0,opts,rp,op,fp,strp);
%     sol.classification = my_classifier(odes_sol.y(:,end));
%   y = odes_sol.y(:,end);
%   
%   keyboard
%
%   tic;
%   kAprev = 2*fp.kA;
%   kA = kAprev/2;
%   y0 = [Bmax,1,1,1];
%   maxits = 100;
%   numits = 0;
%   while (abs(kA-kAprev)/kA > tol) && (numits < maxits)
%     numits = numits+1;
%     [y0,fvals,flag,output] = ...
%       fsolve(@(y) dydt_for_eq(0,y,rp,op,new_fp(kA),strp),y0,...
%       optimoptions('fsolve','MaxFunctionEvaluations',1e5,...
%       'MaxIterations',1e5,'Display','off'));
%     kAprev = kA;
%     kA = blendfact*(1/(1-percentage)-1)/y0(3)+(1-blendfact)*kAprev;
%     kAs(numits) = kA
%   end
%   toc
% %   
% %   tic;
% %   %  y0 = [Bmax,1,1,1,fp.kA];
%   y0 = [Bmax,1,1,1,6.667];
%   y = fsolve(@(y) dydt_for_eq2(0,y,rp,op,new_fp(y(5)),strp),y0,...
%     optimoptions('fsolve','MaxFunctionEvaluations',1e5,...
%     'MaxIterations',1e5,'SpecifyObjectiveGradient',true,...
%     'FunctionTolerance',1e-15,'OptimalityTolerance',1e-15,...
%     'Display','off'));
% %   toc
% 
%   tic;
  maxits = 1e4;
  success = false;
  ycurrs = [Bmax 0;1 1;Bmax 0];
  ytmp = [];
  for eqc = 1:2
    ycurr = ycurrs(:,eqc);
    yprev = ycurr*2;
    numits = 0;
    shrink_factor = sqrt(tol^(1/maxits));
    while (max(abs(ycurr-yprev)./ycurr) > tol) && (numits < maxits)
      numits = numits+1;
      [f,J] = dydt_for_eq3(0,ycurr,rp,op,fp,strp);
      yprev = ycurr;
      ys(:,numits) = ycurr;
      dy = J\f;
%       dy = dy.*min(shrink_factor^numits,1);
      delta = min(1,0.9*min(yprev./abs(dy)));
      ycurr = yprev-delta*dy;
      ytmp = [ytmp,ycurr];
    end
    y(:,eqc) = ycurr;
    if (numits >= maxits)
%       ycurr = ycurrs(:,eqc);
      [y(:,eqc),fval(eqc)] = fminsearch(@(y) sum(...
        dydt_for_eq3(0,y,rp,op,fp,strp).^2),...
        ycurr,optimset('MaxFunEvals',1e5,'MaxIter',1e5,...
        'TolX',1e-15,'TolFun',1e-15));
    else
      fval(eqc) = sum(f.^2);
    end
  end
  if (y(1) < 1)
    if fval(1) < fval(2)
      ycurr = y(:,1);
    else
      ycurr = y(:,2);
    end
  else
    ycurr = y(:,1);
    if fval(1) > tol
      warning('Issues with finding kA!\n');
    end
  end
  
  y = ycurr;
  frac = 1-percentage;  %  Should match with dydt_for_eq2
  y(4) = (fp.f*y(2)-op.T)*frac*op.tau;
  y(3) = (rp.kAMEps*(rp.a1*frac+rp.sA)*y(4)+...
    (y(2)+fp.nu4)*rp.sA+y(2)*rp.a1*frac)/...
    ((y(4)*rp.kAMEps+y(2)+fp.nu4)*rp.muA);
  y(5) = ycurr(3);
  y(6) = (1-frac)/y(3)/frac;
%   toc
  
%   warning('Solved for kA!');
  fp.kA = y(end);
  fp.Bmax = y(1);
%   fprintf('kA = %30.20g\n',fp.kA);
%   fp = setfield(fp,'kA',y(5));
%   ind = strmatch('kA',klabels,'exact');
%   k(ind) = y(5);

end
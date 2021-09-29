function [values,isterminal,direction] = myevent(t,y,rp,op,fp,strp,...
  tolvec,tstart)
  
  values(1,1) = any(abs(dydt(t,y,rp,op,fp,strp)) > ...
    max(tolvec,tolvec./eps.*eps(y'))');%*ones(size(y));
  values(1,1) = values(1,1) || ...
    (toc(tstart) < min([op.minintegrationtime]));
  
  values(2,1) = toc(tstart) < max([op.maxintegrationtime]);
  
  isterminal = true(size(values)); %true(size(y));
  
  direction = zeros(size(values));%zeros(size(y));

end
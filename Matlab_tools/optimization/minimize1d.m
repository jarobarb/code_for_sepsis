function [graphs] = minimize1d(fnct,df,a,b,c1,c2,option,messup)

tol = 100*sqrt(eps);
[atemp,btemp,c] = findbds(fnct,a,b,1);
if strcmp(c,'Too Big')
    [atemp,btemp,c] = findbds(fnct,a,b,-1);
    if strcmp(c,'Too Big')
        error('Can''t bound any minimums or maximums, look at boundaries.');
        return;
    elseif (feval(fnct,a)-feval(fnct,b))*(b-a)>0
        [atemp,btemp,c] = findbds(fnct,atemp,btemp);
    else
        [atemp,btemp,c] = findbds(fnct,btemp,c);    
        if strcmp(c,'Too Big')
            error('Can''t find minimum');
            return;
        end
    end
end
% atemp,btemp,c

if messup == 1
    messup = .5-rand;
else
    messup = 0;
end

if strcmp(option,'Golden Section')
    [graphs] = goldsec(atemp,btemp+messup*(max(atemp,c)-btemp),c,fnct,tol);
elseif strcmp(option,'Brent')
    [graphs] = brent(atemp,btemp+messup*(max(atemp,c)-btemp),c,fnct,tol);
elseif strcmp(option,'DBrent')
    [graphs] = dbrent(min(atemp,c),btemp+messup*(max(atemp,c)-btemp),max(atemp,c),fnct,df,tol);
elseif strcmp(option,'GenAlg')
    [graphs] = genalg1d(c1,c2,fnct,tol);
elseif strcmp(option,'simann1d')
    [graphs] = simann1d(a,b,c2-c1,fnct,tol);
end
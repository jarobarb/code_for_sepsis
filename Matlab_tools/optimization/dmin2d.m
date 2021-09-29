function [p,fp,counttotal] = min2d(p,f,df,tol)

fp = feval(f,p);
n = length(p);
xi = feval(df,p);
maxits = 200;
counttotal = 0;

g = xi;
h = g;
xi = h;

for its = 1:maxits
    [a0,b0,c0] = findbdsdir(f,0,1,xi,p,1);
    [mini,fmin,count] = brentdir(a0,b0,c0,f,df,xi,p,tol);
    counttotal = counttotal+count;
    fret = fmin;
    xi = mini*xi;
    p = xi+p;
    if 2*abs(fret-fp)<=tol*(abs(fret)+abs(fp)+eps)
        return
    end
    fp = feval(f,p);
    xi = feval(df,p);
    gg = sum(g.^2);
%     dgg = sum(xi.^2);
    dgg = sum((xi+g).*xi);
    if norm(xi) <= eps
        return
    end
    gam = dgg/gg;
    g = -xi;
    h = g+gam*h;
    xi = h;
    fp;
end
error('Maxits exceeded');
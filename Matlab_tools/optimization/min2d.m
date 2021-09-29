function [graphs,p,fp,counttotal,fail] = min2d(p,f,df,tol)

fp = feval(f,p);
n = length(p);
xi = feval(df,p);
maxits = 200;
counttotal = 0;
graphs = [counttotal;fp;p];
small = 1;

g = -xi;
h = g;
xi = h;

for its = 1:maxits
    [a0,b0,c0,fail] = findbdsdir(f,0,small,xi,p,1);
    if fail == 1
        return;
    end
    [mini,fmin,count] = brentdir(a0,b0,c0,f,xi,p,tol);
    counttotal = counttotal+count;
    fret = fmin;
    xi = mini*xi;
    p = xi+p;
    graphs  = [graphs,[counttotal;fmin;p]];
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
fail = 1;
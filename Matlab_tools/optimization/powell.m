function [graphs,p,fp,counttotal] = powell(p,xi,f,tol)

fret = feval(f,p);
pt = p;
n = length(p);
small = 1;
iter = 0;
counttotal = 0;
done = 0;
maxits = 200;
graphs = [0;fret;p];

while ~(done)
    iter    =   iter+1;
    fp      =   fret;
    ibig    =   0;
    del     =   0;
    for ii = 1:n
        xit     =   xi(:,ii);
        fptt    =   fret;
        [a0,b0,c0] = findbdsdir(f,0,small,xit,p,1);
        [mini,fmin,count] = brentdir(a0,b0,c0,f,xit,p,tol);
        fret = fmin;
        xit = mini*xit;
        p = xit+p;
        counttotal = counttotal+count;
        graphs = [graphs,[counttotal;fmin;p]];
        if abs(fptt-fret) > del
            del = abs(fptt-fret);
            ibig=ii;
        end
    end
    if 2*abs(fp-fret) <= tol*(abs(fp)+abs(fret))
        return;
    elseif iter == maxits
        error('Maxits exceeded.');
        return;
    end
    ptt = 2*p-pt;
    xit = p-pt;
    pt = p;
    fptt = feval(f,ptt);
    if fptt < fp
        t = 2*(fp-2*fret+fptt)*(fp-fret-del)^2-del*(fp-fptt)^2;
        if t < 0
            [a0,b0,c0] = findbdsdir(f,0,small,xit,p,1);
            [mini,fmin,count] = brentdir(a0,b0,c0,f,xit,p,tol);
            fret = fmin;
            xit = mini*xit;
            p = xit+p;
            graphs = [graphs,[counttotal;fmin;p]];
            if mini>tol
                xi(:,ibig) = xit*(1/mini);
            end
        end
    end
end
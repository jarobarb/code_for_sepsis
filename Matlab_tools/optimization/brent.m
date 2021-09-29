function [graph] = brent(a0,b0,c0,fnct,tol)

%  Initialize frequently used variables.
maxits = 10;
C = (3-sqrt(5))/2;
zerobuffer = sqrt(eps);
count = 0;
graph = zeros(3,0);

% global colorcount
% global symbolcount
% colors = 'rb';
% symbols = 'x+';

%  Initialize a, b, and c by makin a<b<c (apparently we ought to do this).
a = min(a0,c0);
b = max(a0,c0);
v = b0;
w = v;
x = v;
disttrav = 0;

%  Initialize some functional values.
fx = feval(fnct,x);
fv = fx;
fw = fx;

%  Initially store some test variables.
xm = .5*(a+b);
tol1 = tol*abs(x)+zerobuffer;
tol2 = 2.*tol1;

%  Do the loop while we haven't looped too many times or we've found the
%  answer.
while (count <= maxits)&&(abs(x-xm) > tol2-.5*(b-a))
    
    %  If we're still far from the solution, try a parabolic fit:
    if abs(disttrav)>tol1
        r=(x-w)*(fx-fv);
        q=(x-v)*(fx-fw);
        p=(x-v)*q-(x-w)*r;
        q=2*(q-r);
        if q > 0
            p=-p;
        end
        q=abs(q);
        disttravtemp=disttrav;
        disttrav=d;
        
        %  So we've stored some of the parabolic information.  Now we test
        %  its various qualities.  Hmm, 
        if (abs(p)>=abs(.5*q*disttravtemp))||(p<=q*(a-x))||(p>=q*(b-x))
            [d,disttrav] = block1(x,xm,disttrav,a,b,C);
            [u,fu,a,b,v,fv,w,fw,x,fx,count]=block2(d,tol1,x,fx,w,fw,a,b,v,fv,fnct,count);
        else
            d = p/q;
            u = x+d;
            if (u-a<tol2)||(b-u<tol2)
                d = tol1*sign(xm-x);
            end
            [u,fu,a,b,v,fv,w,fw,x,fx,count]=block2(d,tol1,x,fx,w,fw,a,b,v,fv,fnct,count);
        end
    else
        [d,disttrav] = block1(x,xm,disttrav,a,b,C);
        [u,fu,a,b,v,fv,w,fw,x,fx,count]=block2(d,tol1,x,fx,w,fw,a,b,v,fv,fnct,count);
    end    
    xm = .5*(a+b);
    tol1 = tol*abs(x)+zerobuffer;
    tol2 = 2.*tol1;
    graph = [graph,[count;x;fx]];
%     hold on
%     semilogy(count,abs(fx-1),'o');
end
mini=x;
fmin=fx;
if (count>maxits)
    error('Maximum number of iterations exceeded.');
end

function [d,e] = block1(x,xm,e,a,b,C)
if x>=xm
    e = a-x;
else
    e = b-x;
end
d=C*e;

function [u,fu,a,b,v,fv,w,fw,x,fx,count]=block2(d,tol1,x,fx,w,fw,a,b,v,fv,fnct,count)
if abs(d)>=tol1
    u = x+d;
else
    u = x+tol1*sign(d);
end
fu=feval(fnct,u);
count = count + 1;
if fu <= fx
    if u >= x
        a = x;
    else
        b=x;
    end
    v=w;
    fv=fw;
    w=x;
    fw=fx;
    x=u;
    fx=fu;
else
    if (u < x)
        a = u;
    else
        b = u;
    end
    if (fu <= fw)||(w==x)
        v = w;
        fv = fw;
        w = u;
        fw = fu;
    elseif (fu<=fv)||(v==x)||(v==w)
        v = u;
        fv = fu;
    end
end
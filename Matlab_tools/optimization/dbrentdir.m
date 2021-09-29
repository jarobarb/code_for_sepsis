function [mini,fmin,count] = dbrentdir(a0,b0,c0,f,df,dirr,spot,tol);

maxits=100;  zerobuffer=10*eps;  count=0;  exit=0;

a = min(a0,c0);
b = max(a0,c0);
v = [b0,feval(f,b0*dirr+spot),feval(df,b0*dirr+spot)'*dirr];
w = v;
x = v;
u = v;
e = 0;

xm = .5*(a+b);
tol1=tol*abs(x(1))+zerobuffer;
tol2=2.*tol1;

while ((count < maxits) && (abs(x(1)-xm)>tol2-.5*(b-a))) && (exit == 0)
    if abs(e)>tol1
        d1 = 2*(b-a);
        d2 = d1;
        if (w(3)~=x(3))
            d1=(w(1)-x(1))*x(3)/(x(3)-w(3));
        end
        if v(3)~=x(3)
            d2=(v(1)-x(1))*x(3)/(x(3)-v(3));
        end
        u1=x(1)+d1;
        u2=x(1)+d2;
        test1=((a-u1)*(u1-b)>0) && (x(3)*d1<=0);
        test2=((a-u2)*(u2-b)>0) && (x(3)*d2<=0);
        olde = e;
        e = d;
        if test1 && test2
            if abs(d1)<abs(d2)
                d=d1;
            else
                d=d2;
            end
        elseif test1
            d=d1;
        elseif test2
            d=d2;
        else
            [e,d] = block1(x,a,b);
        end
        if (abs(d)>abs(.5*olde)) && ~(test1 || test2)
            [e,d] = block1(x,a,b);
        end
        if ~(test1 || test2)
            u(1)=x(1)+d;
            if (u(1)-a<tol2) || (b-u(1)<tol2)
                d=sign(xm-x(1))*tol1;
            end
        end
        [u,a,b,v,w,x,exit,count] = block2(d,tol1,x,a,b,w,v,f,df,exit,count,dirr,spot);
    else
        [e,d] = block1(x,a,b);
        [u,a,b,v,w,x,exit,count] = block2(d,tol1,x,a,b,w,v,f,df,exit,count,dirr,spot);
    end
    xm = .5*(a+b);
    tol1=tol*abs(x(1))+zerobuffer;
    tol2=2.*tol1;
end
if count>=maxits
    display('dbrent had too many iterations!');
end
mini=x(1);
fmin=x(2);

function [e,d] = block1(x,a,b)
if x(3)>=0
    e=a-x(1);
else
    e=b-x(1);
end
d=.5*e;

function [u,a,b,v,w,x,exit,count] = block2(d,tol1,x,a,b,w,v,f,df,exit,count,dirr,spot)
if abs(d) >= tol1
    u(1)=x(1)+d;
    u(2)=feval(f,u(1)*dirr+spot);
else
    u(1)=x(1)+tol1*sign(d);
    u(2)=feval(f,u(1)*dirr+spot);
    if u(2)>x(2)
        exit=1;
    end
end
u(3)=feval(df,u(1)*dirr+spot)'*dirr;
count = count+2;
if (u(2)<=x(2)) && (exit == 0)
    if (u(1)>x(1)) || (u(1)==x(1) && floor(count/4)==count/4)
        a=x(1);
    else
        b=x(1);
    end
    v=w;
    w=x;
    x=u;
else
    if u(1)<x(1)
        a=u(1);
    else
        b=u(1);
    end
    if (u(2)<=w(2)) || (w(1)==x(1))
        v=w;
        w=u;
    elseif (u(2)<=v(2)) || (v(1) == x(1)) || (v(1) == w(1))
        v=u;
    end
end
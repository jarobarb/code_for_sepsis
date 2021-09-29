function [graph] = goldsec(a,b,c,fnct,tol)

%  We note that tol should be of sqrt(epsilon).

%  Define the numbers to be used to split the segments.
C = (3-sqrt(5))/2;
R = 1-C;
count = 0;
graph = zeros(3,0);

%  Block down the interval.
x0 = a;
x3 = c;
if abs(c-b) > abs(b-a)
    x1 = b;
    x2 = b+C*(c-b);
else
    x2 = b;
    x1 = b-C*(b-a);
end
f1=feval(fnct,x1);
f2=feval(fnct,x2);

%  It should be noted that this test is used because (as they note on page
%  391) because sqrt(2*|f(b)|/(b^2*f''(b))) is usually of order unity.
while abs(x3-x0)>tol*(abs(x1)+abs(x2))
    if f2<f1
        x0=x1;
        x1=x2;
        x2=R*x1+C*x3;
        f1=f2;
        f2=feval(fnct,x2);
    else
        x3=x2;
        x2=x1;
        x1=R*x2+C*x0;
        f2=f1;
        f1=feval(fnct,x1);
    end
%     hold on
%     semilogy(count,abs(min(f1,f2)-1),'>');
    count = count+1;
    if (f1<f2)
        graph = [graph,[count;x1;f1]];
    else
        graph = [graph,[count;x2;f2]];
    end
end

%  Store the result.
if (f1<f2)
    fmin=f1;
    mini=x1;
else
    fmin=f2;
    mini=x2;
end
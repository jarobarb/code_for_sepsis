function [A,B,X0] = minlogittopandbottom(xst,ystopt,ysbottomt)

global xs ystop ysbottom;
xs = xst; ystop = ystopt; ysbottom = ysbottomt;
tol = 1e-6;
ABX0vertices = [0,1,0;...
                1,-1,.1;...
                -1,1,.1;...
                0,-1,-.1];

[graphs,vertices,y,iter,fail] = simplex(ABX0vertices,@f,tol)

A = mean(vertices(:,1));
B = mean(vertices(:,2));
X0 = mean(vertices(:,3));

function norm2 = f(rowvec);

A = rowvec(1);
B = rowvec(2);
X0 = rowvec(3);

global xs ystop ysbottom;

analytic = xs;
for j=1:length(xs)
    analytic(j) = myphipsi(xs(j),A,B,X0);
end
plot(xs,analytic,'b',xs,ystop,'r+-',xs,ysbottom,'k+-');
A,B,X0
pause

norm2 = sqrt(sum((ysbottom-analytic).^2))+sqrt(sum((ystop-analytic).^2));
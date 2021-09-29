function [Ab,At,x0,w] = minswapna(xst,yst)

global xs ys;
xs = xst; ys = yst;
tol = 1e-6;
bump = .1;
Ab = 7.2; At = 6.9; x0 = 6; w = 2;
Swapnavertices = [Ab,At,x0,w;...
                Ab+bump,At,x0,w;...
                Ab,At+bump,x0,w;...
                Ab,At,x0+bump,w;...
                Ab,At,x0,w+bump];

[graphs,vertices,y,iter,fail] = simplex(Swapnavertices,@f,tol)

Ab = mean(vertices(:,1));
At = mean(vertices(:,2));
x0 = mean(vertices(:,3));
w = mean(vertices(:,4));

function norm2 = f(rowvec);

Ab = rowvec(1);
At = rowvec(2);
x0 = rowvec(3);
w = rowvec(4);

global xs ys;

analytic = ys;
for j=1:length(xs)
    analytic(j) = mylorenz(xs(j),Ab,At,x0,w);
end
plot(xs,analytic,'b',xs,ys,'r');
Ab,At,x0,w
pause

norm2 = sqrt(sum((ys-analytic).^2));
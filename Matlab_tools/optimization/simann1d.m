function [graphs] = simann1d(a0,c0,bounds,f,tol)

count = 1;
bounds = bounds/2;
r = round(rand);
b0 = r*a0+(1-r)*c0;
maxits = 9001;
l = logspace(0,floor(log10(tol)),maxits);
graphs = zeros(3,0);

sol = b0;
currentf = feval(f,sol);
mini = sol;
fmin = currentf;
colors = 'rbk';

global colorcount

while count <= maxits
    newsol = sol+(randn*bounds);
    newf = feval(f,newsol);
    count = count + 1;
    if (newf < currentf) ||...
            (newf > currentf && rand < exp(((currentf-newf)/(abs(currentf+newf))*floor(count/2))))
        sol = newsol;
        currentf = newf;
        if currentf < fmin
            fmin = currentf;
            mini = sol;
        end
    end
    if mod(count,30)-1 == 0
        graphs = [graphs,[count;mini;fmin]];
    end
end
function [graphs,mini,fmin,count] = simann2d(a0,c0,f,tol,scale)

count = 1;
maxits = 9001;
l = logspace(0,floor(log10(tol)),maxits);
graphs = zeros(2,0);

sol = [a0;c0];
currentf = feval(f,sol);
mini = sol;
fmin = currentf;
colors = 'rbk';

global colorcount

while count <= maxits
    direction = exp(i*rand*2*pi);
    direction = randn*scale*[real(direction);imag(direction)];
    newsol = sol+direction;
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
    if mod(count,30)-1==0
        graphs = [graphs,[count;fmin]];
    end
%     figure(3);
%     semilogy(count,abs(fmin),['o',colors(colorcount)])
    hold on
end
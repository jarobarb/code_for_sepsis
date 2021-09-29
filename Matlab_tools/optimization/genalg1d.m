function [graphs] = genalg1d(a0,c0,f,tol)

gensize = 30;    %  Size of Generation.
discret = ceil(-log2(tol));    %  Size of discretization.
fitness = zeros(gensize,1);     %  Stores fitness.
reprodu = zeros(gensize,1);     %  Stores reproduction "wheel" of fortune.
generat = zeros(gensize,discret);   %  Stores generation representatives.
realval = zeros(gensize,1);     %  Stores actual values of binary repres.
preprod = 1;
pmutat = .05;                    %  These two are important variables, perhaps.
mutcha = .7;                   %  Different values for these guys would be good.
count = 1;
maxits = gensize*discret*15;
mini = a0;
fmin = feval(f,a0);
bump = 1;
colors = 'rbk';
graphs = zeros(3,0);

global colorcount;

for members = 1:gensize
    generat(members,:) = round(rand(1,discret));     %  Store 1st gen.
end
generation = 1;

while (count <= maxits)
    
    %  Calculate the floating point values of these binary expansions:
    for members = 1:gensize
        realval(members) = (sum(2.^(discret-1:-1:0).*generat(members,:)))/(2^discret-1);
    end
    realval = a0+realval.*(c0-a0);
    
    %  Calculate the respective function values and shift them properly in
    %  order for a greater or lesser death rate.
    for members = 1:gensize
        fitness(members) = feval(f,realval(members));
        count = count + 1;
    end
        
    %  Store functional values.
    if fmin > min(fitness)
        [tf,loc] = ismember(min(fitness),fitness);
        mini = realval(loc);
        fmin = fitness(loc);
    end
    
    %  Make reproduction wheel.
    fitnessnew = zeros(size(fitness));
    for members = 1:gensize
        if ~(ismember(fitness(members),fitness(1:members-1)))
            fitnessnew = fitnessnew + (fitness < fitness(members));
        end
    end
    fitness = fitnessnew + bump;
    
    %  Normalize reproduction wheel.
    fitness = fitness./(sum(union(fitness,fitness(1))));
    ordered = union(fitness,fitness(1));
    reprodu(1) = fitness(1);
    for members = 2:gensize
        if ismember(fitness(members),fitness(1:members-1))
            reprodu(members) = reprodu(members-1);
        else
            reprodu(members) = fitness(members)+reprodu(members-1);
        end
    end
       
    %  Define genepool.
    pool = rand(gensize,1);
    newpool = pool;
    for members = 1:gensize
        [tf,loc] = ismember(pool(members),union(reprodu,pool(members)));
        newpool(members) = loc;
    end
    pool = newpool;
    
    %  Simple check to make sure we can pair all members of our genepool.
    if gensize/2 ~= round(gensize/2)
        error('Generation size not even!');
        return;
    end
    
    %  Get the next generation of breeders.
    oldgenerat = generat;
    for members=1:gensize
        oldgenerat(members,:) = generat(pool(members),:);
    end
    generat = oldgenerat;
    
    %  Crossover (sometimes depending on preprod).
    for pairs = 1:gensize/2
        if rand < preprod
            crossnum = round(rand*(discret-2))+2;
            generat(2*pairs-1,crossnum:discret) = oldgenerat(2*pairs,crossnum:discret);
            generat(2*pairs,crossnum:discret) = oldgenerat(2*pairs-1,crossnum:discret);
        end
    end

%     for pairs = 1:gensize/2
%         if rand < preprod
%             for members = 1:discret
%                 if rand<.5
%                     generat(2*pairs-1,members) = oldgenerat(2*pairs,members);
%                     generat(2*pairs,members) = oldgenerat(2*pairs-1,members);
%                 else
%                     generat(2*pairs-1,members) = oldgenerat(2*pairs-1,members);
%                     generat(2*pairs,members) = oldgenerat(2*pairs,members);
%                 end
%             end
%         end
%     end
    
    %  Mutate them.
    for members = 1:gensize
        if rand < pmutat
            generat(members,:) = mod(generat(members,:)+round(mutcha*rand(1,discret)),2);
        end
    end
        
%     semilogy(generation,abs(fmin-1),['x',colors(colorcount)])
%     hold on
    graphs = [graphs,[count;mini;fmin]];
    generation = generation + 1;
end
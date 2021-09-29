function [allgraphs,graphtotals,table1,table2] = test2doptimization()

methods = {'Simplex';...
        'Direction Set';...
        'Conjugate Gradient';...
        'SimAnn';...
        'GenAlg'};
functions = {@para;@dpara;@revpeak;@drevpeak};
exactfmin = [1,-8.10621358944099];
exacttol = [1.1,-6.10621358944099];
table = cell(6,4);
table{1,2} = 'count'; table{1,3} = 'fminave'; table{1,4} = 'coef'; table{1,5} = 'successes';
successes = zeros(2,5);
graphlimit = zeros(2,3);
graphcrit = 1;
noofruns = 100;
trialrun = 10;
adj = [1 1 1 1 1];

for j = 2:6
    table{j,1} = methods{j-1};
end
table1 = table;
table2 = table;
graphtotals = cell(2,5);
allgraphs = cell(2,5,noofruns);
fendaverage = zeros(2,5);
counts = zeros(2,5);
counts2 = cell(2,3);
counts3 = cell(2,3);
rand('state',sum(100*clock));
r = [4*rand(noofruns,1)-2,4*rand(noofruns,1)-2]
pause

colsym1 = {'b<','r>';'mx','k+';'bs','ro';'ms','ko';'b<','r>'};

for method = 1:length(methods)
    for fncts = 1:2
        for runs = 1:noofruns
            if method == 1
                [graphs,a,b,c,fail] = simplex([r(runs,:);[4*rand-2,4*rand-2];[4*rand-2,4*rand-2]],...
                    functions{fncts*2-1},10^-6);
            elseif method == 2
                [graphs,a,b,c] = powell(r(runs,:)',eye(2),functions{fncts*2-1},10^-6);
            elseif method == 3
                [graphs,a,b,c,fail] = min2d(r(runs,:)',functions{fncts*2-1},functions{fncts*2},10^-6);
            elseif method == 4
                [graphs] = simann2d(r(runs,1),r(runs,2),functions{fncts*2-1},10^-6,4);
            else
                [graphs,xm,fm,co] = genalg2d([-2;2],[-2;2],functions{fncts*2-1},10^-6);
            end
            [method,fncts,runs]
            [m,n] = size(graphs);
            if (mod(runs,trialrun)-1) == 0
                if runs == 1
                    tempf = 0;
                else
                    tempf = tempf+fendaverage(fncts,method);
                end
                fendaverage(fncts,method)=graphs(2,n);
            else
                fendaverage(fncts,method)=min(graphs(2,n),fendaverage(fncts,method));
                if runs == noofruns
                    tempf = tempf+fendaverage(fncts,method);
                end
            end
            if graphs(2,n) < exacttol(fncts)
                successes(fncts,method) = successes(fncts,method)+1;
                n2 = length(graphs(2,:)) - length(graphtotals{fncts,method});
                allgraphs{fncts,method,runs} = graphs;
                if method <= 3
                    if length(graphtotals{fncts,method}) == 0
                        counts2{fncts,method} = 1:length(graphs(2,:));
                        counts3{fncts,method} = length(graphs(2,:));
                        errors = [abs(graphs(2,1:length(graphs(2,:))-1)-exactfmin(fncts));...
                                abs(graphs(2,2:length(graphs(2,:)))-exactfmin(fncts))];
                    else
                        counts2{fncts,method} = [counts2{fncts,method},1:length(graphs(2,:))];
                        counts3{fncts,method} = [counts3{fncts,method},length(graphs(2,:))];
                        errors = [errors, [abs(graphs(2,1:length(graphs(2,:))-1)-exactfmin(fncts));...
                                abs(graphs(2,2:length(graphs(2,:)))-exactfmin(fncts))]];
                    end                        
                    graphtotals{fncts,method} = [graphtotals{fncts,method},graphs(:,:)];
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                    if successes(fncts,method) == graphcrit
                        graphlimit(fncts,method) = length(graphtotals{fncts,method}(1,:));
                    end
                elseif length(graphtotals{fncts,method}) == 0
                    graphtotals{fncts,method} = graphs;
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                elseif mod(runs,trialrun) ~= 0
                    graphtotals{fncts,method}(2,:) = min([graphtotals{fncts,method}(2,:);...
                            graphs(2,:)],[],1);
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                elseif runs == trialrun
                    temp = graphtotals{fncts,method}(2,:);
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                else
                    temp = graphtotals{fncts,method}(2,:)+temp;
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                end
            end
        end
        counts(fncts,method) = counts(fncts,method)/successes(fncts,method);
        if method > 3
            graphtotals{fncts,method}(2,:) = temp./(noofruns/trialrun);
        end
        [m,n] = size(graphtotals{fncts,method});
        if fncts == 1
            table1{method+1,2} = counts(fncts,method);
            table1{method+1,3} = tempf/(noofruns/trialrun);
            table1{method+1,5} = successes(fncts,method);
        else
            table2{method+1,2} = counts(fncts,method);
            table2{method+1,3} = tempf/(noofruns/trialrun);
            table2{method+1,5} = successes(fncts,method);
        end
        if method == 1
            figure(1);
            hold on
            semilogy(graphtotals{fncts,method}(1,1:graphlimit(fncts,method)),...
                abs(graphtotals{fncts,method}(2,1:graphlimit(fncts,method))-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        elseif method <= 3
            figure(2);
            hold on
            semilogy(graphtotals{fncts,method}(1,1:graphlimit(fncts,method)),...
                abs(graphtotals{fncts,method}(2,1:graphlimit(fncts,method))-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        else
            figure(3);
            hold on
            semilogy(graphtotals{fncts,method}(1,:)./adj(method),...
                abs(graphtotals{fncts,method}(2,:)-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        end
        if method <= 3
            p = polyfit(counts2{fncts,method}(1,:),...
                log(abs(graphtotals{fncts,method}(2,:)-exactfmin(fncts))+eps),1);
        else
            p = polyfit(graphtotals{fncts,method}(1,:)./adj(method),...
                log(abs(graphtotals{fncts,method}(2,:)-exactfmin(fncts))+eps),1);
        end
        if fncts == 1
            table1{method+1,4} = p(1);
        else
            table2{method+1,4} = p(1);
        end            
    end
end

figure(1);
legend('Simplex-Fnct1','Simplex-Fnct2');
figure(2);
legend('Direction-Fnct1','Direction-Fnct2','Conjugate-Fnct1','Conjugate-Fnct2');
figure(3);
legend('SimAnn-Fnct1','SimAnn-Fnct2','GenAlg-Fnct1','GenAlg-Fnct2');
format long g
table1
table2
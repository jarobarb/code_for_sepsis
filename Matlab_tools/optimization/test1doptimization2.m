function [allgraphs,graphtotals,table1,table2] = test1doptimization2()

methods = {'Golden Section';...
        'Brent';...
        'DBrent';...
        'simann1d';...
        'GenAlg'};
functions = {@fnct1d1;
    @dfnct1d1;
    @fnct1d2;
    @dfnct1d2};
exactfmin = [1;-96.50140856037186294];
exacttol = [1.1,-50];
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
r = [rand(1,noofruns);rand(1,noofruns)]
pause
b1 = 0; b2 = 1;

colsym1 = {'b<','r>';'mx','k+';'bs','ro';'ms','ko';'b<','r>'};

for method = 1:length(methods)
    for fncts = 1:2
        for runs = 1:noofruns
            [method fncts runs]
            [graphs] = minimize1d(functions{fncts*2-1},...
                functions{fncts*2},r(1,runs),r(2,runs),b1,b2,methods{method},fncts);
            [m,n] = size(graphs);
            if (mod(runs,trialrun)-1) == 0
                if runs == 1
                    tempf = 0;
                else
                    tempf = tempf+fendaverage(fncts,method);
                end
                fendaverage(fncts,method) = graphs(3,n);
            else
                fendaverage(fncts,method) = min(fendaverage(fncts,method),graphs(3,n));
            end
            if runs == 100
                tempf = tempf+fendaverage(fncts,method);
            end
            if graphs(3,n) < exacttol(fncts)
                successes(fncts,method) = successes(fncts,method)+1;
                n2 = length(graphs(3,:)) - length(graphtotals{fncts,method});
                allgraphs{fncts,method,runs} = graphs;
                if method <= 3
                    if length(graphtotals{fncts,method}) == 0
                        counts2{fncts,method} = 1:length(graphs(3,:));
                        counts3{fncts,method} = length(graphs(3,:));
                        errors = [abs(graphs(3,1:length(graphs(3,:))-1)-exactfmin(fncts));...
                                abs(graphs(3,2:length(graphs(3,:)))-exactfmin(fncts))];
                    else
                        counts2{fncts,method} = [counts2{fncts,method},1:length(graphs(3,:))];
                        counts3{fncts,method} = [counts3{fncts,method},length(graphs(3,:))];
                        errors = [errors, [abs(graphs(3,1:length(graphs(3,:))-1)-exactfmin(fncts));...
                                abs(graphs(3,2:length(graphs(3,:)))-exactfmin(fncts))]];
                    end                        
                    graphtotals{fncts,method} = [graphtotals{fncts,method},graphs(:,:)];
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                    if successes(fncts,method) == graphcrit
                        graphlimit(fncts,method) = length(graphtotals{fncts,method}(1,:));
                    end
                elseif length(graphtotals{fncts,method}) == 0
                    graphtotals{fncts,method} = graphs;
%                     graphtotals{fncts,method}(3,:) =...
%                         log(abs(graphtotals{fncts,method}(3,:)-exactfmin(fncts))+eps);
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                elseif mod(runs,trialrun) ~= 0
                    graphtotals{fncts,method}(3,:) = min([graphtotals{fncts,method}(3,:);...
                            graphs(3,:)],[],1);
%                     graphtotals{fncts,method}(3,:) = log(abs(graphs(3,:)-exactfmin(fncts))+eps)+...
%                         graphtotals{fncts,method}(3,:)-graphs(3,:);
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                elseif runs == trialrun
                    temp = graphtotals{fncts,method}(3,:);
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                else
                    temp = graphtotals{fncts,method}(3,:)+temp;
                    counts(fncts,method) = counts(fncts,method)+graphs(1,n);
                end
            end
        end
        counts(fncts,method) = counts(fncts,method)/successes(fncts,method);
        if method > 3
            graphtotals{fncts,method}(3,:) = temp./(noofruns/trialrun);
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
            figure(4);
            hold on
            semilogy(graphtotals{fncts,method}(1,1:graphlimit(fncts,method)),...
                abs(graphtotals{fncts,method}(3,1:graphlimit(fncts,method))-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        elseif method <= 3 && method ~= 1
            figure(4);
            hold on
            semilogy(graphtotals{fncts,method}(1,1:graphlimit(fncts,method)),...
                abs(graphtotals{fncts,method}(3,1:graphlimit(fncts,method))-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        else
            figure(5);
            hold on
            semilogy(graphtotals{fncts,method}(1,:)./adj(method),...
                abs(graphtotals{fncts,method}(3,:)-exactfmin(fncts))+eps,...
                colsym1{method,fncts});
        end
        if method <= 3
            p = polyfit(counts2{fncts,method}(1,:),...
                log(abs(graphtotals{fncts,method}(3,:)-exactfmin(fncts))+eps),1);
        else
            p = polyfit(graphtotals{fncts,method}(1,:)./adj(method),...
                log(abs(graphtotals{fncts,method}(3,:)-exactfmin(fncts))+eps),1);
        end
        if fncts == 1
            table1{method+1,4} = p(1);
        else
            table2{method+1,4} = p(1);
        end            
    end
end

figure(4);
legend('Golden Section-Fnct1','Golden Section-Fnct2','Brent-Fnct1','Brent-Fnct2','DBrent-Fnct1','DBrent-Fnct2');
figure(5);
legend('SimAnn-Fnct1','SimAnn-Fnct2','GenAlg-Fnct1','GenAlg-Fnct2');
format long g
table1
table2
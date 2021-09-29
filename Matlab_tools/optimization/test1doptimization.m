function [graphtotals,table1,table2] = test1doptimization()

methods = {'Golden Section';...
        'Brent';...
        'DBrent';...
        'GenAlg';...
        'simann1d'};
functions = {@fnct1d1;
    @dfnct1d1;
    @fnct1d2;
    @dfnct1d2};
actualfmin = [1;-96.50140856037186294];
table = cell(6,5);
table{1,2} = 'mini'; table{1,3} = 'fmin'; table{1,4} = 'count'; table{1,5} = 'coef';
R = rand(1,10);

for j = 2:6
    table{j,1} = methods{j-1};
end 
table1 = table;
table2 = table;

colsym1 = {'b<:','r>:';'bx-.','r+-.';'bs--','ro--';'bs--','ro:';'b<--','r>:'};

for method = 1:length(methods)
    for fncts = 1:2
        if method <= 3
            figure(4);
            [graphs] = minimize1d(functions{fncts*2-1},...
                functions{fncts*2},.5-rand,1,methods{method});
            [m,n] = size(graphs);
            if fncts == 1
                table1{method+1,2} = graphs(2,n);
                table1{method+1,3} = graphs(3,n);
                table1{method+1,4} = graphs(1,n);
            else
                table2{method+1,2} = graphs(2,n);
                table2{method+1,3} = graphs(3,n);
                table2{method+1,4} = graphs(1,n);
            end 
            if fncts == 1
                hold on
                semilogy(graphs(1,:),abs(graphs(3,:)-actualfmin(1))+eps,colsym1{method,fncts});
                p = polyfit(graphs(1,:),log(abs(graphs(3,:)-actualfmin(1))+eps),1);
                table1{method+1,5} = p(1);
            else
                hold on
                semilogy(graphs(1,:),abs(graphs(2,:)-actualfmin(2))+eps,colsym1{method,fncts});
                p = polyfit(graphs(1,:),log(abs(graphs(2,:)-actualfmin(2))+eps),1);                
                table2{method+1,5} = p(1);
            end 
        else
            figure(5);
            for ii = 1:10
%                 colorcount = ii;
                if ii == 1
                    [graphs] = minimize1d(functions{fncts*2-1},...
                        functions{fncts*2},-1,2,methods{method});
                else
                    graphs = graphs+minimize1d(functions{fncts*2-1},...
                        functions{fncts*2},-1,2,methods{method});
                end 
            end 
            graphs = graphs.*(.1);
            [m,n] = size(graphs);
            if fncts == 1
                table1{method+1,2} = graphs(2,n);
                table1{method+1,3} = graphs(3,n);
                table1{method+1,4} = graphs(1,n);
            else
                table2{method+1,2} = graphs(2,n);
                table2{method+1,3} = graphs(3,n);
                table2{method+1,4} = graphs(1,n);
            end 
            if fncts == 1
                hold on
                semilogy(graphs(1,:),abs(graphs(3,:)-actualfmin(1)),colsym1{method,fncts});
                p = polyfit(graphs(1,:),log(abs(graphs(3,:)-actualfmin(1))+eps),1);
                table1{method+1,5} = p(1);
            else
                hold on
                semilogy(graphs(1,:),abs(graphs(3,:)-actualfmin(2)),colsym1{method,fncts});
                p = polyfit(graphs(1,:),log(abs(graphs(3,:)-actualfmin(2))+eps),1);                
                table2{method+1,5} = p(1);
            end 
        end             
    end 
end 

figure(4);
legend('GoldSec-Fnct1','GoldSec-Fnct2','Brent-Fnct1','Brent-Fnct2','DBrent-Fnct1','DBrent-Fnct2');
figure(5);
legend('Genalg-Fnct1','Genalg-Fnct2','SimAnn-Fnct1','SimAnn-Fnct2');
table1
table2

% for ii = 1:3
%     for method = 4:5
%         [mini,fmin,count] = minimize1d(inline('(.535*x-.23).^2'),...
%             inline('2*(.535*x-.23)'),0,1,methods{method});
%     end
% end
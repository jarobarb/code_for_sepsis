function FRBC = myphipsi(FB,A,B,X0)
    FRBC = FB;
    for j=1:length(FB)
        if (FB(j) <= X0)
            FRBC(j) = 0;
        elseif (FB(j) >= 1-X0)
            FRBC(j) = 1;
        else
            Q = (FB(j)-X0)/(1-2*X0);
            FRBC(j) = 1/(exp(-A)*(Q/(1-Q))^(-B)+1);
        end
    end
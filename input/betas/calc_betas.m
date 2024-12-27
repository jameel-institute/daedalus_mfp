addpath('../');

diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Delta','Covid Wildtype','SARS'};
beta_fixed = zeros(size(diseases));

for j = 1:length(diseases);    

    inp2 = diseases(j);

    if strcmp(inp2,'Influenza 2009');
        dis = param_influenza_2009;
    elseif strcmp(inp2,'Influenza 1957');
        dis = param_influenza_1957;
    elseif strcmp(inp2,'Influenza 1918');
        dis = param_influenza_1918;
    elseif strcmp(inp2,'Covid Omicron');
        dis = param_covid_omicron;    
    elseif strcmp(inp2,'Covid Delta');
        dis = param_covid_delta;
    elseif strcmp(inp2,'Covid Wildtype');
        dis = param_covid_wildtype;    
    elseif strcmp(inp2,'SARS');
        dis = param_sars;
    else
        error('Unknown Disease!');
    end  

    T1 = readtable(strcat('LLMIC_',string(inp2),'.csv'));
    T2 = readtable(strcat('UMIC_',string(inp2),'.csv'));
    T3 = readtable(strcat('HIC_',string(inp2),'.csv'));
    T  = [T1;T2;T3]; 

    r0a           = T.res1;
    betabar       = dis.R0/mean(r0a);%this is the beta such that mean(betabar*r0a) = R0
    beta_fixed(j) = betabar;  

    % betas   = T.res2;
    % betabar = mean(betas);
    % disp(betabar*mean(r0a)-dis.R0);%not necessarily the same!

end
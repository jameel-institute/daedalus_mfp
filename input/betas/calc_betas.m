addpath('../');

diseases = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Wildtype','Covid Omicron','Covid Delta','SARS'};
beta_fix = zeros(size(diseases));

for j = 1:length(diseases);    

    inp2 = diseases(j);

    if strcmp(inp2,'Influenza 2009');
        dis = param_influenza_2009;
    elseif strcmp(inp2,'Influenza 1957');
        dis = param_influenza_1957;
    elseif strcmp(inp2,'Influenza 1918');
        dis = param_influenza_1918;
    elseif strcmp(inp2,'Covid Wildtype');
        dis = param_covid_wildtype;    
    elseif strcmp(inp2,'Covid Omicron');
        dis = param_covid_omicron;    
    elseif strcmp(inp2,'Covid Delta');
        dis = param_covid_delta;    
    elseif strcmp(inp2,'SARS');
        dis = param_sars;
    else
        error('Unknown Disease!');
    end  

    T1 = readtable(strcat('HIC_',string(inp2),'.csv'));
    T2 = readtable(strcat('UMIC_',string(inp2),'.csv'));
    T3 = readtable(strcat('LLMIC_',string(inp2),'.csv'));
    T  = [T1;T2;T3]; 

    R01         = T.GDPL;
    betastar    = dis.R0/mean(R01);%this is the beta such that mean(betastar*R01) = R0
    beta_fix(j) = betastar;  

    %betas   = T.VSYL;
    %betabar = mean(betas);
    %disp(betabar*mean(R01)-dis.R0);%not necessarily the same!

end
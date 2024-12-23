function sec = p2Sim(inp1,inp2,inp3)
    
    load('countries/country.mat','data');
    CD        = readtable('countries/country_data.csv');
    data      = p2RandIG(data,CD,inp1);
    data.tvec = 1+[0 365*4];
    
    [dis,p2] = p2Params(data,inp2);
    
    data = p2Strat(data,inp3);

    [data,f,g] = p2Run(data,dis,p2);

    [cost,~,~] = p2Cost(data,dis,g);
    sec        = sum(cost([3,6,7:end],:),'all');
    %sec(1)    = sum(cost([3],:),'all');
    %sec(2)    = sum(cost([7:end],:),'all');
    %sec(3)    = sum(cost([6],:),'all');

    p2Plot(data,f,p2,g,cost,NaN,inp1,inp2,inp3);
    
end
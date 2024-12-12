function sec = p2Sim(inp1,inp2,inp3)
    
    load(strcat('countries/',inp1,'.mat'),'data');
    data.tvec = 1+[0 365*3];
    
    %[data,~,~]    = p2Params(data,'Covid Wildtype');%to define wnorm and Td_CWT
    [data,dis,p2] = p2Params(data,inp2);
    
    data = p2Strat(data,inp3);

    [data,f,g] = p2Run(data,dis,inp3,p2);

    [cost,~] = p2Cost(data,dis,p2,g);
    sec(1)   = sum(cost([3,6,7:10],:),'all');
    sec(2)   = sum(cost([3],:),'all');
    sec(3)   = sum(cost([6],:),'all');
    sec(4)   = sum(cost([7:10],:),'all');

    p2Plot(data,f,p2,g,cost,NaN,sec(1),inp1,inp2,inp3);
    
end
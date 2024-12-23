function sec = dd_single(inp1,inp2,inp3)
    
    load('input/country.mat','data');
    CD        = readtable('input/country_data.csv');
    data      = dd_set_country(data,CD,inp1);
    data.tvec = 1+[0 365*4];
    
    [dis,p2] = dd_set_disease(data,inp2);
    
    data = dd_set_strategy(data,inp3);

    [data,f,g] = dd_run_sim(data,dis,p2);

    [cost,~,~] = dd_calc_loss(data,dis,g);
    sec        = sum(cost([3,6,7:end],:),'all');
    %sec(1)    = sum(cost([3],:),'all');
    %sec(2)    = sum(cost([7:end],:),'all');
    %sec(3)    = sum(cost([6],:),'all');

    dd_plot_output(data,f,p2,g,cost,NaN,inp1,inp2,inp3);
    
end
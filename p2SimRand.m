function p2SimRand
    
    parpool;
    
    taskdir = strcat('output/',datestr(now,'yyyy-mm-dd_HH-MM-SS'),'/');%'./';
    mkdir(taskdir);
    
    locations  = {'LLMIC','UMIC','HIC'};
    diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Delta','Covid Wildtype','SARS'};
    strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
    
    lloc       = length(locations);
    nsamples   = 5000;
    samples    = table();
    ldis       = length(diseases);
    lstrat     = length(strategies);
    
    data_array = cell(lloc,nsamples,lstrat);
    dis_array  = cell(lloc,nsamples,ldis);
    p2_array   = cell(lloc,nsamples,ldis);

    load('countries/country.mat','data');
    data.tvec    = 1+[0 365*4];
    country_data = readtable('countries/country_data.csv');
    
    for h = 1:lloc;
        inp1 = locations{h};
        parfor i = 1:nsamples;
            ldata   = data;
            ldata   = p2RandIG(ldata,country_data,inp1);
            %for l  = 1:2;if l==2;ldata.t_vax = min(ldata.t_vax,100);end
            row     = p2InpTable(inp1,i,ldata);
            samples = [samples;row];
            for j = 1:ldis;
                inp2             = diseases{j};    
                [dis,p2]         = p2Params(ldata,inp2);
                dis_array{h,i,j} = dis;
                p2_array{h,i,j}  = p2;
            end
            for k = 1:lstrat;
                inp3              = strategies{k};
                ldata             = p2Strat(ldata,inp3);
                data_array{h,i,k} = ldata;
            end
        end
    end
    writetable(samples, strcat(taskdir,'sample_data.csv'));
    
    for h = 1:lloc;
    for j = 1:ldis;
    for k = 1:lstrat;
        output = [];
        parfor i = 1:nsamples;
        data = data_array{h,i,k};
        dis  = dis_array{h,i,j};
        p2   = p2_array{h,i,j};
        try
            [~,~,g] = p2Run(data,dis,p2);
            [~,~,c] = p2Cost(data,dis,g);
            sec     = c;
        catch
            sec      = nan(1,16);
        end
        output = [output;sec];
        end
        p2OutTable(output,strcat(string(locations{h}),'_',string(diseases{j}),'_',string(strategies{k})),taskdir);
        %p2Plot(data,f,p2,g,cost,NaN,inp1,inp2,inp3);
        %disp([h,j,k,i]);
    end
    end
    end
    
    delete(gcp);
    
end
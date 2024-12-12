function p2SimRand%(locations,diseases,strategies)
    
    parpool;
    taskdir = './';%strcat('output/',datestr(now,'yyyy-mm-dd_HH-MM-SS'),'/');%task directory
    %mkdir(taskdir);

    locations  = {'United States'};
    diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Wildtype','Covid Omicron','Covid Delta','SARS'};
    strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
    
    llocs     = length(locations);
    nsamples  = 1000;%samples per location
    ldis      = length(diseases);
    lstrat    = length(strategies);
    inputs    = zeros(llocs,nsamples,493);
    dat_array = cell(llocs,nsamples,ldis,lstrat);
    dis_array = cell(llocs,nsamples,ldis,lstrat);
    p2_array  = cell(llocs,nsamples,ldis,lstrat);
    CD        = readtable('countries/country_data.csv');
    load('countries/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    data.tvec = 1+[0 365*3];
    
    for h = 1:llocs;
        
        inp1 = locations{h};
                
        parfor i = 1:nsamples;
    
        ldata         = data;
        ldata         = p2RandCountry(ldata,CD,inp1);
        %ldata        = p2RandIG(ldata,CD,inp1);
        inputs(h,i,:) = [i,NaN,NaN,ldata.Npop',ldata.NNs(1:45)',...
                         ldata.CM(:)',ldata.comm,ldata.travelA3,ldata.schoolA1,ldata.schoolA2,ldata.workp,...
                         ldata.obj',ldata.wfh(1,:),ldata.wfh(2,:),...
                         ldata.t_vax,ldata.arate,ldata.puptake,ldata.Hmax,ldata.t_tit,ldata.trate,ldata.Tres,ldata.sdl,ldata.sdb,...
                         NaN,ldata.la];
        #[ldata,~,~]   = p2Params(ldata,'Covid Wildtype');%to define wnorm and Td_CWT
        %for l = 1:2;if l==2;ldata.t_vax = min(ldata.t_vax,100);end
                   
        for j = 1:ldis;
        
        inp2           = diseases{j};    
        [ldata,dis,p2] = p2Params(ldata,inp2);
        
        for k = 1:lstrat;
        
        inp3  = strategies{k};
        ldata = p2Strat(ldata,inp3);
    
        dat_array{h,i,j,k} = ldata;
        dis_array{h,i,j,k} = dis;
        p2_array{h,i,j,k}  = p2;
    
        end
        end
        end
    end
    
    for j = 1:ldis;
        for k = 1:lstrat;
        for h = 1:llocs;

        outputs = zeros(nsamples,4);

        parfor i = 1:nsamples;
        
        ldata = dat_array{h,i,j,k};
        dis   = dis_array{h,i,j,k};
        inp3  = strategies{k};
        p2    = p2_array{h,i,j,k};

        try
            [~,~,g]  = p2Run(ldata,dis,inp3,p2);
    
            [cost,~] = p2Cost(ldata,dis,p2,g);
            sec(1)   = sum(cost([3,6,7:10],:),'all');
            sec(2)   = sum(cost([3],:),'all');
            sec(3)   = sum(cost([6],:),'all');
            sec(4)   = sum(cost([7:10],:),'all');
        catch
            sec = nan(1,4);
        end
        outputs(i,:) = sec;
        %disp([h,i,j,k]);
        %plots = p2Plot(data,f,p2,g,cost,ccost_t,sec(1),inp1,inp2,inp3);

        end

        p2OutTable([squeeze(inputs(h,:,:)),outputs],strcat(string(locations{h}),'_',string(diseases{j}),'_',string(strategies{k})),taskdir);

        end
        end
    end

    delete(gcp);
    
end
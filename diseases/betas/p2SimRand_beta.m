function p2SimRand_beta(igroups,diseases)%,strategies)
    
    parpool;
    taskdir = strcat('diseases/',datestr(now,'yyyy-mm-dd_HH-MM-SS'),'/');%task directory
    mkdir(taskdir);
    
    ligroups  = length(igroups);
    nsamples  = 3000;%samples per income group
    CD        = readtable('countries/country_data.csv');
    load('countries/Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields    = fieldnames(data);
    ikeep     = [6,7,8,13,14,16,17,18];
    data      = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    data.tvec = 1+[0 365*3];
    ldis      = length(diseases);
    lstrat    = 1;%length(strategies);
    
    for h = 1:ligroups;
    
        inp1    = igroups{h};
        inputs  = zeros(nsamples,493);
        outputs = zeros(ldis,lstrat,nsamples,4);
        
        parfor i = 1:nsamples;
    
        ldata       = data;
        ldata       = p2RandCountry(ldata,CD,inp1);
        inputs(i,:) = [i,NaN,NaN,ldata.Npop',ldata.NNs(1:45)',...
                       ldata.CM(:)',ldata.comm,ldata.travelA3,ldata.schoolA1,ldata.schoolA2,ldata.workp,...
                       ldata.obj',ldata.wfh(1,:),ldata.wfh(2,:),...
                       ldata.t_vax,ldata.arate,ldata.puptake,ldata.Hmax,ldata.t_tit,ldata.trate,ldata.Tres,ldata.sdl,ldata.sdb,...
                       NaN,ldata.la];
        [ldata,~,~] = p2Params(ldata,'Covid Wildtype');%to define wnorm and Td_CWT
        %for l = 1:2;if l==2;ldata.t_vax = min(ldata.t_vax,100);end
                   
        for j = 1:ldis;
        
        inp2           = diseases{j};    
        [ldata,dis,p2] = p2Params(ldata,inp2);
        
        for k = 1:lstrat;
        
%         inp3  = strategies{k};
%         ldata = p2Strat(ldata,inp3);
%         
%         try
%             [~,~,g]  = p2Run(ldata,dis,inp3,p2);
%     
%             [cost,~] = p2Cost(ldata,dis,p2,g);
%             sec(1)   = sum(cost([3,6,7:10],:),'all');
%             sec(2)   = sum(cost([3],:),'all');
%             sec(3)   = sum(cost([6],:),'all');
%             sec(4)   = sum(cost([7:10],:),'all');
%         catch
%             disp([h,i,j,k]);
%             sec = nan(1,4);
%         end
% 
%         outputs(j,k,i,:) = sec;
%         %plots = p2Plot(data,f,p2,g,cost,ccost_t,sec(1),inp1,inp2,inp3);
                
        R0a = dis.R0/dis.beta;
        outputs(j,k,i,:) = [NaN,NaN,dis.beta,R0a];     

        end
        end
        end
    
        for j = 1:ldis;
        for k = 1:lstrat;
        p2OutTable([inputs,squeeze(outputs(j,k,:,:))],strcat(string(igroups{h}),'_',string(diseases{j})),taskdir);
        end
        end

    end
    
    delete(gcp);
    
end
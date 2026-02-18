function main_r0rt
    
parpool;

taskdir = './';
%mkdir(taskdir);
addpath('../../');

locations  = {'LLMIC','UMIC','HIC'};
diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918','Covid Omicron','Covid Delta','Covid Wildtype'};
strategies = {'School Closures','School Closures','Economic Closures','Economic Closures','Elimination'};

lloc       = length(locations);
nsamples   = 5000;
ldis       = length(diseases);
lstrat     = length(strategies);

data_array = cell(lloc,nsamples,lstrat);
dis_array  = cell(lloc,nsamples,ldis);
p2_array   = cell(lloc,nsamples,ldis);

load('input/country.mat','data');
data.tvec    = 1+[0 365*10];%
country_data = readtable('input/country_data.csv');

for h = 1:lloc;
    inp1    = locations{h};
    samples = table();
    parfor i = 1:nsamples;
        ldata   = data;
        ldata   = dd_set_country(ldata,country_data,inp1);%for l  = 1:2;if l==2;ldata.t_vax = min(ldata.t_vax,100);end
        row     = dd_store_input(i,ldata);
        samples = [samples;row];
        for j = 1:ldis;
            inp2             = diseases{j};    
            [dis,p2]         = dd_set_disease(ldata,inp2);
            dis_array{h,i,j} = dis;
            p2_array{h,i,j}  = p2;
        end
        for k = 1:lstrat;
            inp3              = strategies{k};
            ldata             = dd_set_strategy(ldata,inp3);
            data_array{h,i,k} = ldata;
        end
    end
    %writetable(samples,strcat(taskdir,string(locations{h}),'_data.csv'));
end

for h = 1:lloc;
for j = 1:ldis;
for k = 1:lstrat;
    output = [];
    parfor i = 1:nsamples;
    data = data_array{h,i,k};
    dis  = dis_array{h,i,j};
    p2   = p2_array{h,i,j};
    %
    rt   = calc_r0rt(data,dis,p2,k);
    %
    output = [output;[dis.r0,rt]];
    end
    T = array2table(output);
    T.Properties.VariableNames = {'r0','rt'};
    writetable(T,strcat(taskdir,strcat(string(locations{h}),'_',string(diseases{j}),'_',string(strategies{k}),'_',string(k)),'.csv'));
    %dd_store_output(output,strcat(string(locations{h}),'_',string(diseases{j}),'_',string(strategies{k}),'_',string(k))),taskdir);
    %p2Plot(data,g,p2,f,cost,inp1,inp2,inp3);
    %disp([h,j,k,i]);
end
end
end

delete(gcp);

end
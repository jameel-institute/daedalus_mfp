%% inputs

parpool;

warning('off');

filename  = '../../../Data/Preparedness/0.income_group.xlsx';
T1        = readtable(filename);
filename  = '../../../Data/Preparedness/1.population_age.xlsx';
T2        = readtable(filename);
countries = intersect(T1.Var1,T2.Location);
data      = struct;
% filename  = '../../../Data/Preparedness/5.gva_sector.xlsx';
% T         = readtable(filename);
% countries = T.ToIndustry_Sector;

%% loop

%for i = 170;%:length(countries);
parfor i = 1:length(countries);

country = countries{i};

%% income group

filename = '../../../Data/Preparedness/0.income_group.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Var1,country));

%year   = T(kr,:).Var2;
igroup = char(T(kr,:).Var3);
%gnipc  = T(kr,:).Var40;

if strcmp(igroup,'LIC')||strcmp(igroup,'LMIC');
igroup = 'LLMIC';
end

%data.igroup = igroup;
%data.gnipc  = gnipc;

%disp([country,' is a ',igroup,' and its GNI per capita is $',num2str(gnipc),' (World Bank, ',num2str(year),')']);

%% population by age

filename = '../../../Data/Preparedness/1.population_age.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Location,country));

% if ~isempty(kr);
kc   = find(strcmp(T.Properties.VariableNames,'x0_4'));
Npop = 1000*table2array(T(kr,kc:kc+20))';
% else
%     continue;
%     % Npop = zeros(21,1);
% end

%data.Npop = Npop;

%disp(['The population is ',num2str(sum(Npop)),' (United Nations, 2019)']);

%% population by sector
  
filename = '../../../Data/Preparedness/2.population_sector.xlsx';
T        = readtable(filename);
kr       = find(strcmp(T.Country_,country));

if ~isempty(kr);
    kc     = find(strcmp(T.Properties.VariableNames,'Agriculture_Hunting_Forestry'));        
    year   = T(kr,:).Year_;
    NNs    = table2array(T(kr,kc:kc+44))';
    NEC    = T(kr,:).NotElsewhereClassified;
    source = 'OECD/ILO';
    %% informal employment    
    filename = '../../../Data/Preparedness/2.population_sector_informal.xlsx';
    T        = readtable(filename);
    kr       = find(strcmp(T.Country,country));
    
    if ~isempty(kr);
        kc   = find(strcmp(T.Properties.VariableNames,'x4_A_Agriculture_ForestryAndFishing'));
        Ninf = table2array(T(kr,kc:kc+20))';
        ind  = [1,1,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,...
                4,5,6,7,8,8,8,8,8,9,10,10,10,11,12,13,14,15,16,17,18,19,20]';
        Ntot = accumarray(ind,NNs);
        disc = max(0,Ninf(1:20)-Ntot);
        for k = 1:45;
            indj   = ind(k);
            if Ntot(indj)>0;
                prop = NNs(k)/Ntot(indj);
            else
                prop = 1/sum(ind==indj);
            end
            NNs(k) = NNs(k) + disc(indj)*prop; 
        end
        NNs = NNs + max(NEC,Ninf(21))*NNs/sum(NNs);       
    else                
        NNs = NNs + NEC*NNs/sum(NNs); 
    end
    %% scaling
    if year<2019;
        filename = '../../../Data/Preparedness/2.population_sector_change.xlsx';
        T        = readtable(filename);
        T.Properties.VariableNames = {'a','b','c','d','x2009','x2010','x2011','x2012','x2013',...
                                      'x2014','x2015','x2016','x2017','x2018','x2019'};
                                  
        kr   = find(strcmp(T.b,country));
        kc   = find(strcmp(T.Properties.VariableNames,['x',num2str(year)]));
        dpop = T(kr,:).x2019/table2array(T(kr,kc));
        NNs  = dpop*NNs;
    end
    %% capping
    NNs = round(NNs);%redistribution of NEC and scaling may both result in non-integer NNs
    if sum(NNs)>sum(Npop(5:13));
        NNs = sum(Npop(5:13))*NNs/sum(NNs);
        NNs = round(NNs);
        while sum(NNs) > sum(Npop(5:13));
            ind         = find(NNs == max(NNs)); 
            NNs(ind(1)) = NNs(ind(1)) - 1; 
        end
    end
else
    NNs = [];%continue;
    % year   = 0;
    % NNs    = zeros(45,1);
    % source = 'Estimated';
end

%NNs = [NNs;Npop(1);sum(Npop(2:4));sum(Npop(5:13))-sum(NNs);sum(Npop(14:end))];

%data.NNs = NNs;

%disp(['The workforce accounts for ',num2str(round(100*sum(NNs(1:45))/sum(Npop(5:13)))),...
%      '% of the adult population, across ',num2str(nnz(NNs(1:45))),' sectors (',source,', ',num2str(year),')']);

%% contact matrices

filename = '../../../Data/Preparedness/3.contact_matrices_all.xlsx';

if any(strcmp(sheetnames(filename),country));
    opts       = detectImportOptions(filename);
    opts.Sheet = country;
    
    filename = '../../../Data/Preparedness/4.contact_rates_home.xlsx';
    T        = readtable(filename,opts);
    matAL    = table2array(T);
    matAL    = compress_mat(matAL,Npop);

    filename = '../../../Data/Preparedness/4.contact_rates_other.xlsx';
    T        = readtable(filename,opts);
    matAHT   = table2array(T);
    matAHT   = compress_mat(matAHT,Npop);

    filename = '../../../Data/Preparedness/4.contact_rates_school.xlsx';
    T        = readtable(filename,opts);
    matAS    = table2array(T);
    matAS    = compress_mat(matAS,Npop);
    
    %if ~isempty(NNs);
    filename = '../../../Data/Preparedness/4.contact_rates_work.xlsx';
    T        = readtable(filename,opts);
    matBC    = table2array(T);
    workp    = dot(sum(matBC(5:13,:),2),Npop(5:13))/sum(Npop(5:13));%average number of workplace contacts per adult (working and non-working)
    %*(sum(Npop(5:13))/sum(NNs(1:45)));%average number of contacts per working adult
    %else
    %workp    = [];
    %end
    
    source = 'Prem et al., 2021';
else
    matAL  = [];%continue;
    matAHT = [];
    matAS  = [];
    workp  = [];

    % source = 'Estimated, 0';
end

% data.matAL  = matAL;
% data.matAHT = matAHT;
% data.matAS  = matAS;
% data.workp  = workp;

%disp(['The average number of contacts per person per day is ',...
%      num2str(round(dot(sum(CM,2),[Npop(1:15);sum(Npop(16:end))])/sum(Npop),2)),' (',source,')']);

%% gva by sector

filename = '../../../Data/Preparedness/5.gva_sector.xlsx';
T        = readtable(filename);

kr = find(strcmp(T.ToIndustry_Sector,country));
if ~isempty(kr); 
    kc = find(strcmp(T.Properties.VariableNames,'D01T02_Agriculture_Hunting_Forestry'));
    
    obj    = abs(table2array(T(kr,kc:kc+44))');%one value negative for some reason
    source = 'OECD, 2018';
    %% scaling
    filename = '../../../Data/Preparedness/5.gva_sector_change.xlsx';
    T        = readtable(filename);
    
    kr = find(strcmp(T.CountryName,country));
    
    dgva = T(kr,:).x2019/T(kr,:).x2018;
    if (dgva>0 && dgva<Inf);
        obj(obj>=0) = dgva*obj(obj>=0);
    end
else
    obj = [];%continue;
    % obj    = zeros(45,1);
    % source = 'Estimated, 0';
end

%data.obj = obj;

%disp(['The GDP is $',num2str(sum(365*obj)),' million (',source,')']);

%% economic closures

% filename   = '../../../Data/Preparedness/6.economic_closures.xlsx';
% opts       = detectImportOptions(filename);
% 
% opts.Sheet = 'Australia 45';
% T          = readtable(filename,opts);
% x_elim     = T.Var2;
% x_elim(41) = 1.00;
% 
% opts.Sheet   = 'United Kingdom 45';
% T            = readtable(filename,opts);
% x_econ(:,1)  = T.Aug;
% x_econ(:,2)  = T.Apr;
% x_econ(41,1) = 1.00;
% x_econ(41,2) = 0.10;
% 
% opts.Sheet   = 'Indonesia 45';
% T            = readtable(filename,opts);
% x_schc(:,1)  = T.TriwulanIV;
% x_schc(:,2)  = T.TriwulanII;
% x_schc(41,1) = 0.10;
% x_schc(41,2) = 0.10;
% 
% % data.alp    = alp;
% %data.x_elim = x_elim;
% %data.x_econ = x_econ;
% %data.x_schc = x_schc;
% 
% %disp(['The labour share of productivity is ',num2str(alp),' (',source,')']);

%% home-working

filename   = '../../../Data/Preparedness/7.home_working.xlsx';
opts       = detectImportOptions(filename);
opts.Sheet = 'LD vs. FO';
T          = readtable(filename,opts);
w_uk       = T.Var5(2:end);

if ~isempty(NNs);
    opts.Sheet = 'Gottlieb';
    T          = readtable(filename,opts);
    kr         = find(strcmp(T.Var1,country));
    if ~isempty(kr);
        w      = T(kr,:).Var2;
    else
        opts.Sheet         = 'GDPpcppp';
        opts.VariableNames = arrayfun(@num2str,1:68,'UniformOutput',0);
        opts.DataRange     = 'A4:BP270';
        T                  = readtable(filename,opts);
        kr                 = find(strcmp(T.x1,country));
        gdppcppp           = str2double(cell2mat(T(kr,:).x68));%available for everywhere but Taiwan, which doesn't have NNs data either
        w                  = max(0,-0.2444+0.0791*log10(gdppcppp));
    end
    scal    = fzero(@(scal) dot(min(scal*w_uk, 1),NNs(1:45)) - w*sum(NNs(1:45)), 1);%initial guess for scal of 1
    scal    = max(0,scal);%this prevents scal of negative epsilon
    wfh     = min(scal*w_uk, 1)';
    wfh(41) = 0;%treating education sector separately
    source  = 'Gottlieb et al., 2021';
else 
    wfh = [];
end

%data.wfh = wfh;

%disp(['The share of home-working is ',num2str(w),' (',source,')']);

%% vaccination

filename = '../../../Data/Preparedness/8.vaccination.csv';
T        = readtable(filename);

kr = strcmp(T.location,country);
if any(kr);
    %independent variable
    date = T(kr,:).date;
    dmat = datevec(date);
    yvec = dmat(:,1)-2020;
    dvec = day(date,'dayofyear') + 366*min(1,yvec) + 365*max(0,yvec-1);%days since 1st Jan 2020
    %dependent variables
    v1   = T(kr,:).people_vaccinated_per_hundred;%cumulative first doses per 100
    v2   = T(kr,:).people_fully_vaccinated_per_hundred;%cumulative second doses per 100
    % Tts  = T(kr,:).people_vaccinated;
    % Tts  = max(fillmissing(Tts,'linear'),0);
    % v2   = T(kr,:).people_fully_vaccinated;
    % v2   = max(fillmissing(v2,'linear'),0);
    % Vts  = (Tts+v2)/2/sum(Npop/10^5);
    %missing values and combination
    time   = [min(dvec):1:max(dvec)]';
    V1     = nan(size(time));    
    V2     = nan(size(time));
    [~,iv] = ismember(dvec, time);
    V1(iv) = v1;
    V2(iv) = v2;
    V1     = max(0, fillmissing(V1,'linear'));
    V2     = max(0, fillmissing(V2,'linear'));
    Vts    = 1000*(V1+V2)/2;%average cumulative doses per 100k%this is for administration time and rate only!!!
    
    rng default;
    m       = (Vts(end)-Vts(find(Vts>0,1)))/(time(end)-time(find(Vts>0,1)));%approximate slope for IC
    lb      = [300, 365, 0];%parameters are t1, t2 and maximum, which are converted to administration time and rate (uptake calculated separately below)
    x0      = [max(300,time(end)-Vts(end)/m), time(end), Vts(end)];
    ub      = [1000, time(end), 100000];
    fun     = @(params,time)pw_function(params,time);
    options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1000000);
    % problem = createOptimProblem('lsqcurvefit','options',options,'x0',x0,'lb',lb,'ub',ub,'objective',fun,'xdata',time,'ydata',Vts);
    % ms      = MultiStart;
    % poptim  = run(ms,problem,1);
    poptim  = lsqcurvefit(fun, x0, time, Vts, lb, ub, options);
    
    t_vax   = poptim(1);
    arate   = poptim(3)/(poptim(2)-poptim(1));
    puptake = max(T(kr,:).people_fully_vaccinated_per_hundred)/100;%uptake based on second dose only
    puptake = puptake/(1-(1/2.87));%normalise by HIT using R0 from https://pmc.ncbi.nlm.nih.gov/articles/PMC7657547/#:~:text=The%20estimated%20summary%20reproductive%20number,CI%2C%202.39%E2%80%933.44). 
    source  = 'Our World in Data, 2022';

    % figure;
    % hold on;
    % plot(time,Vts);
    % plot(time,pw_function(poptim,time));

else
    t_vax   = [];
    arate   = [];
    puptake = [];
    source  = 'Estimated, 0';
end

% data.t_vax   = t_vax;
% data.arate   = arate;
% data.puptake = puptake;

%disp(['The vaccine rollout begins on day ',num2str(t_vax),...
%      ' with an administration rate of ',num2str(arate),' vaccines per 100k per day',...
%      ' and the uptake is ',num2str(100*puptake),'% (',source,')']);

%% hospital capacity

filename = '../../../Data/Preparedness/9.hospital_capacity.xlsx';
opts     = detectImportOptions(filename);

opts.Sheet = 'Beds';
T          = readtable(filename,opts);
kr         = strcmp(T.Var1,country);
if any(kr);
    %year   = T(kr,:).Var2;
    hcap   = 100*T(kr,:).Var3;
    source = 'World Bank/OECD';
else
    %year   = NaN;
    hcap   = [];
    source = 'Estimated';
end

opts.Sheet = 'BOR';
T          = readtable(filename,opts);
kr         = strcmp(T.Var1,country);
if any(kr);
    BOR = T(kr,:).Var3/100;
else
    BOR = 0.85;%assumption
end

% opts.Sheet = 'Covid';
% T          = readtable(filename,opts);
% kr         = strcmp(T.Var1,country);
% if any(kr);
%     Cmax = max(T(kr,:).Var3)/10;
% else
%     Cmax = 0;
% end

Hmax = hcap*(1-BOR);

%data.Hmax  = hcap*(1-BOR);
%data.SHmax = hcap*(1-BOR)*2;

%disp(['The hospital capacity is ',num2str(hcap*(1-BOR)),' beds per 100k (',source,',',num2str(year),')']);

%% testing

filename = '../../../Data/Preparedness/10.testing.csv';
T        = readtable(filename);

kr = strcmp(T.Entity,country);
if any(kr);
    %independent variable
    date = T(kr,:).Date;
    dmat = datevec(date);
    yvec = dmat(:,1)-2020;
    dvec = day(date,'dayofyear') + 366*min(1,yvec) + 365*max(0,yvec-1);%days since 1st Jan 2020
    %dependent variables
    t1   = 100*T(kr,:).CumulativeTotalPerThousand;%cumulative tests per 100k
    %missing values
    time    = [min(dvec):1:max(dvec)]';
    Tts     = nan(size(time));    
    [~,iv]  = ismember(dvec, time);
    Tts(iv) = t1;
    Tts     = max(0, fillmissing(Tts,'linear'));
    
    rng default;
    %bounds and initial condition
    av_daily = movmean(diff(Tts),[0 30]);
    ttit_ub  = find(av_daily>10,1);
    if ~isempty(ttit_ub);
        ttit_ub = time(ttit_ub);
    else
        ttit_ub = time(end);
    end
    lb      = [14, 365, 0];%parameters are t1, t2 and maximum, which are converted to administration time and rate
    x0      = [ttit_ub/2, time(end), Tts(end)];
    ub      = [ttit_ub, time(end), 3500000];%3.5m is max from data
    fun     = @(params,time)pw_function(params,time);
    options = optimoptions(@lsqcurvefit,'MaxFunctionEvaluations',1000000);
    % problem = createOptimProblem('lsqcurvefit','options',options,'x0',x0,'lb',lb,'ub',ub,'objective',fun,'xdata',time,'ydata',Tts);
    % ms      = MultiStart;
    % poptim  = run(ms,problem,1);
    poptim  = lsqcurvefit(fun, x0, time, Tts, lb, ub, options);
    
    t_tit   = poptim(1)/3.48;%normalise by doubling time from https://pmc.ncbi.nlm.nih.gov/articles/PMC7575205/pdf/nihms-1636611.pdf
    trate   = poptim(3)/(poptim(2)-poptim(1));
    source  = 'Our World in Data, 2022';
    
    % figure;
    % hold on;
    % plot(time,Tts);
    % plot(time,pw_function(poptim,time));

    %old way
    % daily = movmean(diff(Tts),[0 30]);
    % i_tit = find(daily>1,1);
    % if ~isempty(i_tit);
    %     i_tit = i_tit+1;
    %     t_tit = fdat(i_tit);
    % else
    %     i_tit = 367;
    %     t_tit = 367; 
    % end
    % t_end = 731;
    % i_end = find(fdat==t_end,1);
    % tspan = t_end-t_tit;
    % Tspan = fTts(i_end)-fTts(i_tit-1);
    % trate = Tspan/tspan;
    
else
    t_tit  = [];
    trate  = [];
    source = 'Estimated, 0';
end

% dur   = 1;
% Ip    = linspace(0,1000,500);
% trate = 110;
% b0    = 2.197;
% b1    = 0.1838;
% b2    = -1.024;
% p3    = (Ip<trate) .*   (1./(1+exp(b0+b1*Ip+b2*log10(trate))))/dur + ...
%         (Ip>=trate).*min(1./(1+exp(b0+b1*Ip+b2*log10(trate))),trate/10^5)/dur;
% p4    = p3;
% plot(Ip,p3);
% hold on;
% plot(trate*ones(1,11),[0:0.1:1],'r');
% plot(Ip,(trate/10^5)*ones(size(Ip)),'r');
% plot(Ip,(Ip.*p3)/trate,'g');%test positivity rate

% data.t_tit = t_tit;
% data.trate = trate;

%disp(['Mass testing begins on day ',num2str(t_tit),...
%      ' with an administration rate of ',num2str(trate),' tests per 100k per day','% (',source,')']);

%% response time & social distancing

filename = '../../../Data/Preparedness/11.response.csv';
T        = readtable(filename);

kc = find(strcmpi(T.Properties.VariableNames,strrep(strrep(strrep(country,' ',''),'-','_'),'''','_')));
if any(kc); 
    date = T.country_name;
    dmat = datevec(date);
    yvec = dmat(:,1)-2020;
    dvec = day(date,'dayofyear') + 366*min(1,yvec) + 365*max(0,yvec-1);%days since 1st Jan 2020

    Bts    = table2array(T(:,kc));
    i_r    = find(Bts>=20,1);%if strategy is no closures, this makes no difference, and if this threshold isn't met we assume Covid closures weren't imposed 
    Tres   = dvec(i_r)/3.48;%normalise by doubling time from https://pmc.ncbi.nlm.nih.gov/articles/PMC7575205/pdf/nihms-1636611.pdf
    source = 'Blavatnik, 2022';
else
    Tres   = [];
    source = 'Estimated, 0';
end

% data.Tres = Tres;

%disp(['The response time is on the ',num2str(Tres),' day of the first year (',source,')']);

filename = '../../../Data/Preparedness/12.mobility.csv';
T        = readtable(filename);
kr       = strcmp(T.Entity,country);

if any(kr) & ~isempty(Tres);
    date = T(kr,:).Day;
    dmat = datevec(date);
    yvec = dmat(:,1)-2020;
    t2   = day(date,'dayofyear') + 366*min(1,yvec) + 365*max(0,yvec-1);%days since 1st Jan 2020
    Rts  = min(T(kr,:).retail_and_recreation,0);
    Gts  = min(T(kr,:).grocery_and_pharmacy,0);
    Tts  = min(T(kr,:).transit_stations,0);
    Mts  = 1+((Rts+Gts+Tts)/3/100);
    
    filename = '../../../Data/Preparedness/12.excess_deaths.csv';
    T        = readtable(filename);
    kr       = strcmp(T.location_name,country);
    exd      = max(1,T(kr,:).mean_value);
    
    filename = '../../../Data/Preparedness/12.deaths.csv';%daily deaths per million (OWID file no longer supported)
    T        = readtable(filename);
    kc       = find(strcmpi(T.Properties.VariableNames,strrep(strrep(strrep(country,' ',''),'-','_'),'''','_')));
    date     = T.date;
    dmat     = datevec(date);
    yvec     = dmat(:,1)-2020;
    t1       = day(date,'dayofyear') + 366*min(1,yvec) + 365*max(0,yvec-1);%days since 1st Jan 2020
    Dts      = table2array(T(:,kc))/10;%daily deaths per 100k
    %Dts     = movmean(table2array(T(:,kc))/10,7);%7-day average
    t1       = t1(~isnan(Dts));
    Dts      = exd*Dts(~isnan(Dts));
    
    [t,i1,i2] = intersect(t1,t2);
    Dts       = Dts(i1);
    Mts       = Mts(i2);
    t         = t - Tres;%days since response time
    
    %customfit3d   = fittype('l + (1-l)/(1+exp(b*log10(d)*exp(-c*t)))',...
    %                        'dependent',{'sd'},'independent',{'d','t'},'coefficients',{'l','b','c'});
    %log10-logistic functions are used elsewhere (CIT), but this function has a value of l + (1-l)/2 at t = 0 and d = 1, without introducing another parameter for the intercept
    customfit3d   = fittype('l + (1-l)*exp(-b*exp(-c*t)*d)',...
                            'dependent',{'sd'},'independent',{'d','t'},'coefficients',{'l','b','c'});%parameters are minimum, and death and time-steepness
    options       = fitoptions(customfit3d);
    options.lower = [0.1, 0, 0];%lower bound of 0.1 since Mts never fell below 0.1 for any country at any time
    options.upper = [min(Mts), 10000, -log(0.5)/7];%upper bound are values achieved during Covid; c half-life cannot be less than 7 days 
    [mdl,gof]     = fit([Dts,t],Mts,customfit3d,options);
    
    % figure;
    % hold on;
    % scatter3(Dts, t, Mts, 'ko');
    % Drange = linspace(0,5,1000);
    % Trange = linspace(0,900,1000);
    % [D,T]  = meshgrid(Drange,Trange);
    % SD     = reshape(feval(mdl,D(:),T(:)),size(D));
    % surf(D,T,SD);
    % zlim([0 1]);
    % colormap(flipud(cool));
    % shading interp;
    % alpha 0.8;
    % view(0,0);
    
    % low     = 1;
    % upp     = 366;
    % lu      = min(Mts);   
    % tit_fun = @(l,b,x) (l-b)+(1-l+b)*(1+((l-1)/(1-l+b))).^(x./10);
    % tit_fit = fit([Dts(t>low&t<upp)],Mts(t>low&t<upp),tit_fun,...
    %                'StartPoint',[0.25,0.00000001],'Lower',[0.10,0],'Upper',[min(lu,0.4),10],...
    %                'Robust','LAR','MaxIter',10*10^3,'MaxFunEvals',10*10^3);
    % figure;
    % hold on;
    % scatter(Dts(t>low&t<upp),Mts(t>low&t<upp),'go');
    % D = linspace(0,10,1000);      
    % plot(D,tit_fun(tit_fit.l,tit_fit.b,D),'b-');
    % xlim([0 10]);
    % ylim([0 1]);
    
    sdl    = mdl.l;
    sdb    = mdl.b;
    sdc    = mdl.c;
    source = 'Our World in Data/Wang et al., 2022';
else
    %continue;
    sdl    = [];
    sdb    = [];
    sdc    = [];
    source = 'Estimated, 0';
end

% data.sdl = sdl;
% data.sdb = sdb;
% data.sdc = sdc;

%disp(['The transmission modifier asymptote is ',num2str(sdl),'(',source,')']);

%% costing - vlyl

filename   = '../../../Data/Preparedness/13.life_expectancy.csv';
T          = readtable(filename);
kr         = strcmp(T.Location,country);

% %value of statistical life
% filename = '../../../Data/Preparedness/13.vly.csv';
% T        = readtable(filename);
% kr       = strcmp(T.DataSource,country);

if any(kr);
    for k = 1:18;
        label  = strcat('AGE',num2str(5*(k-1)),'-',num2str(5*(k-1)+4));
        index  = kr & strcmp(T.Dim2ValueCode,label);
        la(k)  = mean(T(index,:).Value,1);%assuming equal proportion by sex
    end
    % year     = T(kr,:).Var3;
    % gnipcppp = T(kr,:).Var69;
    % vsl      = 100*gnipcppp/10^6;%160*gnipcppp/10^6;%in millions
    % na         = [Npop(1:17)',sum(Npop(18:end))];%length is 18 to match life table
    % lg         = [dot(la(1),na(1))/sum(na(1)),...
    %               dot(la(2:4),na(2:4))/sum(na(2:4)),...
    %               dot(la(5:13),na(5:13))/sum(na(5:13)),...
    %               dot(la(14:end),na(14:end))/sum(na(14:end))];
    % for k = 1:length(lg); 
    %     lgh(k) = sum(1./((1+0.03).^[1:lg(k)]));
    % end    
    % %value of life-year
    % avage    = round(dot(Npop,[2:5:102]')/sum(Npop));
    % label    = strcat('AGE',num2str(avage-mod(avage,5)),'-',num2str(avage-mod(avage,5)+4));  
    % krv      = strcmp(T.Location,country)&strcmp(T.Dim2ValueCode,label);    
    % rlifex   = mean(T(krv,:).Value,1);
    % vly      = vsl/rlifex;
    % vly      = vsl/(dot(lgh,[Npop(1);sum(Npop(2:4));sum(Npop(5:13));sum(Npop(14:end))])/sum(Npop));
    source = 'WHO, 2019';
else
    la     = [];
    %vly    = [];
    source = 'Estimated, 0';
end

%data.la  = la;
%data.vly = vly;
%disp(['The value of a life year is $',num2str(vly),' million (',source,')']);

%% costing - vsyl

% % filename = '../../Data/15.vsy.xlsx';
% % T        = readtable(filename);
% % kr       = strcmp(T.CountryName,country);
% % 
% % if any(kr);
% %     year   = cell2mat(T(kr,:).Year);
% %     gdp    = T(kr,:).MOSTRECENT;
% %     vsy    = 2.02*gdp/sum(Npop(2:4))/10^6;%in millions
% %     source = strcat('World Bank,',num2str(year));
% % else
% %     vsy    = 0;
% %     source = 'Estimated, 0';
% % end
% 
% gdp = 365*sum(obj);
% vsy = 0.5454*gdp/sum(Npop(2:4));%2.02*gdp/sum(Npop(2:4));
% 
% if vsy~=0;
%     source = 'OECD, 2018';
% else
%     source = 'Estimated, 0';
% end
% 
% data.vsy = vsy;
% 
% % %disp(['The value of a school year is $',num2str(vsy),' million (',source,')']);

%% storage

data(i).igroup  = igroup;
data(i).country = country;
data(i).Npop    = Npop';
data(i).NNs     = NNs';
data(i).matAL   = matAL(:)';
data(i).matAHT  = matAHT(:)';
data(i).matAS   = matAS(:)';
data(i).workp   = workp;
data(i).obj     = obj';
data(i).wfh     = wfh;
data(i).t_vax   = t_vax;   
data(i).arate   = arate;   
data(i).puptake = puptake; 
data(i).Hmax    = Hmax;
data(i).t_tit   = t_tit; 
data(i).trate   = trate; 
data(i).Tres    = Tres;
data(i).sdl     = sdl; 
data(i).sdb     = sdb; 
data(i).sdc     = sdc;
data(i).la      = la;  

% cdata = data(j);
% save(strcat(country,'.mat'),'cdata');

if mod(i,10)==0;
    display(i);
end

% if any(isnan([Npop;NNs;CM(:);comm;hospA2;hospA3;hospA4;travelA3;schoolA1;schoolA2;workp;B;C;obj])); 
%     error(['NaN in data for country: ',country]); 
% end
% save(strcat(country,'.mat'),'data');

end

writetable(struct2table(data), 'country_data.csv');

delete(gcp);

%% postprocessing & plotting

% for i = 1:size(candidates,1)
%     
%     hicc(i)  = strcmp(candidates{i,2},'HIC');
%     umicc(i) = strcmp(candidates{i,2},'UMIC');
%     lmicc(i) = strcmp(candidates{i,2},'LMIC');
%     licc(i)  = strcmp(candidates{i,2},'LIC');
%     
% end
% 
% hics  = candidates(hicc,:);
% umics = candidates(umicc,:);
% lmics = candidates(lmicc,:);
% lics  = candidates(licc,:);
% 
% [~,ihic]  = sort([hics{:,3}], 'ascend','MissingPlacement','first'); 
% [~,iumic] = sort([umics{:,3}],'ascend','MissingPlacement','first'); 
% [~,ilmic] = sort([lmics{:,3}],'ascend','MissingPlacement','first'); 
% [~,ilic]  = sort([lics{:,3}], 'ascend','MissingPlacement','first'); 
% 
% hics  = hics(ihic,:);
% umics = umics(iumic,:);
% lmics = lmics(ilmic,:);
% lics  = lics(ilic,:);
% 
% candidates = [lics;lmics;umics;hics];

% names = categorical(candidates(:,1));
% names = reordercats(names,cellstr(candidates(:,1)));
% 
% f  = figure('Units','centimeters','Position',[0 0 45 15]);
% set(f,'defaulttextInterpreter','latex');
% set(f,'defaultAxesTickLabelInterpreter','latex');
% set(f,'defaultLegendInterpreter','latex');
% set(f,'DefaultAxesFontSize',10);
% fs = 12;
% 
% g = gca;
% g.Position = [0.035 0.28 0.93 0.625];
% hold on;
% 
% b = bar(names,cell2mat(candidates(:,3)));
% %plot(names,70*ones(1,length(candidates)),'-','linewidth',1,'color','k');
% 
% %ylim([0 100000]);
% %ylabel('GNI per capita (\$)');
% g.XAxis.TickLength = [0 0];
% xtickangle(45);
% %g.YAxis.Exponent = 3;
% grid on;
% box on;
% set(gca,'FontSize',fs);
% 
% b.FaceColor = 'flat';
% b.CData(1:size(lics,1),:)                                                        = repmat([0 1 1],size(lics,1),1);
% b.CData(size(lics,1)+1:size(lics,1)+size(lmics,1),:)                             = repmat([0 .5 1],size(lmics,1),1);
% b.CData(size(lics,1)+size(lmics,1)+1:size(lics,1)+size(lmics,1)+size(umics,1),:) = repmat([0 0 .5],size(umics,1),1);
% b.CData(size(lics,1)+size(lmics,1)+size(umics,1)+1:end,:)                        = repmat([0.41 0.16 0.38],size(hics,1),1);

%% functions

function mat = compress_mat(mat,pop)

pop(16) = sum(pop(16:end));
pop     = pop(1:16);

mat = [mat(:,1),sum(mat(:,2:4),2),sum(mat(:,5:13),2),sum(mat(:,14:16),2)];%sum of the columns
mat = [mat(1,:);
       pop(2:4)'*mat(2:4,:)/sum(pop(2:4));
       pop(5:13)'*mat(5:13,:)/sum(pop(5:13));
       pop(14:16)'*mat(14:16,:)/sum(pop(14:16))];%weighted average of the rows

end

function S = pw_function(params,dates)
    
    syms t;
    S = piecewise(t<params(1),0,...
                  params(1)<t<params(2),...
                  -params(1)*(params(3)/(params(2)-params(1))) + t*(params(3)/(params(2)-params(1))),...
                  t>params(2),params(3));
                          
    S = subs(S,t,dates);
    S = fillmissing(double(S),'nearest');
    
end
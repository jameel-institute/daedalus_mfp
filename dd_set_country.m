function data = dd_set_country(data,country_data,location)

if any(strcmp(country_data.igroup,location));
    country_data = country_data(strcmp(country_data.igroup,location),:);
elseif any(strcmp(country_data.country,location));
    country_data = country_data(strcmp(country_data.country,location),:);
end

%Npop
nonempind = find(~isnan(country_data.Npop_1));%indices of possible options
randindex = nonempind(randi(length(nonempind)));%index of randomly selected option
randvalue = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'Npop_')};
if any(strcmp(country_data.igroup,location));
    defivalue = 50*10^6*randvalue'/sum(randvalue);
elseif any(strcmp(country_data.country,location));
    defivalue = randvalue';
end
data.Npop = defivalue;

%la
nonempind = find(~isnan(country_data.la_1));
randindex = nonempind(randi(length(nonempind)));
randvalue = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'la_')};
data.la   = randvalue;

%matAL
nonempind  = find(~isnan(country_data.matAL_1));
randindex  = nonempind(randi(length(nonempind)));
randvalue  = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'matAL_')};
defivalue  = reshape(randvalue,4,4);
data.matAL = defivalue;

%matAHT
nonempind   = find(~isnan(country_data.matAHT_1));
randindex   = nonempind(randi(length(nonempind)));
randvalue   = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'matAHT_')};
defivalue   = reshape(randvalue,4,4);
data.matAHT = defivalue;

%matAS
nonempind  = find(~isnan(country_data.matAS_1));
randindex  = nonempind(randi(length(nonempind)));
randvalue  = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'matAS_')};
defivalue  = reshape(randvalue,4,4);
data.matAS = defivalue;

%workp
nonempind  = find(~isnan(country_data.workp));
randindex  = nonempind(randi(length(nonempind)));
randvalue  = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'workp')};
data.workp = randvalue;

%NNs
nonempind = find(~isnan(country_data.NNs_1));
randindex = nonempind(randi(length(nonempind)));
randvalue = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'NNs_')};%number of workers by sector in real country
defivalue = randvalue/sum(country_data{randindex,2+[5:13]});%proportion of adult population by sector in real country
defivalue = defivalue*sum(data.Npop(5:13));%number of workers by sector in artificial country
defivalue = [defivalue,data.Npop(1),sum(data.Npop(2:4)),sum(data.Npop(5:13))-sum(defivalue),sum(data.Npop(14:end))]';
data.NNs  = defivalue;

%obj
nonempind = find(~isnan(country_data.obj_1)&~isnan(country_data.NNs_1));
randindex = nonempind(randi(length(nonempind)));
randvalue = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'obj_')};%daily gva by sector in real country
if any(strcmp(country_data.igroup,location));
    defivalue                   = randvalue./country_data{randindex,startsWith(country_data.Properties.VariableNames, 'NNs_')};%daily gva per worker by sector in real country
    defivalue(isnan(defivalue)) = 0;
    defivalue(isinf(defivalue)) = 0;
    defivalue                   = defivalue'.*data.NNs(1:45);%daily gva by sector in artificial country
elseif any(strcmp(country_data.country,location));
    defivalue                   = randvalue';
end
data.obj  = defivalue;

%wfh
nonempind = find(~isnan(country_data.wfh_1));
randindex = nonempind(randi(length(nonempind)));
randvalue = country_data{randindex,startsWith(country_data.Properties.VariableNames, 'wfh_')};
data.wfh  = randvalue;

%mitigation measures: note that wblrnd inputs are ordered scale and shape, gamrnd inputs are shape and scale
if strcmp(location,'LLMIC');
    data.Tres    = lognrnd(3.0575,0.2935);
    data.sda     = normrnd(0.6621,1.3124);
    data.sdb     = gamrnd(0.3935,1/0.7127);
    data.sdc     = lognrnd(-4.6282,0.9244);
    data.t_tit   = gamrnd(3.1277,1/0.0620);
    data.trate   = lognrnd(2.4965,1.5519);
    data.Hmax    = lognrnd(2.8324,0.8798);
    data.t_vax   = lognrnd(6.3189,0.1295);
    data.arate   = wblrnd(168.5784,1.0946);
    data.puptake = wblrnd(0.5656,1.2555);
elseif strcmp(location,'UMIC');
    data.Tres    = wblrnd(21.4521,8.3960);
    data.sda     = normrnd(1.0320,1.8192);
    data.sdb     = gamrnd(0.3908,1/0.7409);
    data.sdc     = lognrnd(-4.9317,0.7605);
    data.t_tit   = gamrnd(2.8273,1/0.0707);
    data.trate   = lognrnd(4.4726,0.7496);
    data.Hmax    = lognrnd(3.9377,0.7926);
    data.t_vax   = gamrnd(93.3176,1/0.1971);
    data.arate   = gamrnd(3.5626,1/0.0161);
    data.puptake = wblrnd(0.8960,2.5400);
elseif strcmp(location,'HIC');
    data.Tres    = wblrnd(20.2176,8.1795);
    data.sda     = normrnd(0.4492,0.8620);
    data.sdb     = gamrnd(0.3515,1/0.4975);
    data.sdc     = gamrnd(0.5889,1/152.5585);
    data.t_tit   = lognrnd(3.2143,0.5836);
    data.trate   = lognrnd(5.6181,0.8854);
    data.Hmax    = lognrnd(4.2928,0.7379);
    data.t_vax   = lognrnd(6.0571,0.1035);
    data.arate   = wblrnd(383.6611,3.7781);
    data.puptake = wblrnd(1.2097,7.0783);
elseif any(strcmp(country_data.country,location));
    data.Tres    = country_data{1,startsWith(country_data.Properties.VariableNames, 'Tres')};
    data.sda     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sda')};
    data.sdb     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sdb')};
    data.sdc     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sdc')};
    data.t_tit   = country_data{1,startsWith(country_data.Properties.VariableNames, 't_tit')};
    data.trate   = country_data{1,startsWith(country_data.Properties.VariableNames, 'trate')};
    data.Hmax    = country_data{1,startsWith(country_data.Properties.VariableNames, 'Hmax')};
    data.t_vax   = country_data{1,startsWith(country_data.Properties.VariableNames, 't_vax')};
    data.arate   = country_data{1,startsWith(country_data.Properties.VariableNames, 'arate')};
    data.puptake = country_data{1,startsWith(country_data.Properties.VariableNames, 'puptake')};
end

%vly
%Masterman & Viscusi (2018) method, using US 2019 VSL of $10.9m from 
%https://www.transportation.gov/office-policy/transportation-policy/revised-departmental-guidance-on-valuation-of-a-statistical-life-in-economic-analysis
na        = [data.Npop(1:17)',sum(data.Npop(18:end))];%length is 18 to match life table
la        = data.la;
gdp       = 365*sum(data.obj);%in millions USD
gdppc     = gdp/sum(na);
if all(strcmp(country_data.igroup,'LLMIC')) || (all(strcmp(country_data.igroup,'UMIC')) && gdppc < 0.008809);
    vsl   = 10.9*((0.008809/0.060362)^0.85)*(gdppc/0.008809);
elseif (all(strcmp(country_data.igroup,'UMIC')) && gdppc > 0.008809) || all(strcmp(country_data.igroup,'HIC'));
    vsl   = 10.9*((gdppc/0.060362)^0.85);
end
defivalue = vsl/(dot(la,na)/sum(na));
data.vly  = defivalue;

%vsy
%Psacharopoulos (2021) method, weighting LLMIC value by number of LICs and LMICs, scaled to 1 year
if all(strcmp(country_data.igroup,'LLMIC'));
    llpc  = dot([0.62,0.22],[27,55])/(27+55)/0.33;
elseif all(strcmp(country_data.igroup,'UMIC'));
    llpc  = 0.22/0.33;
elseif all(strcmp(country_data.igroup,'HIC'));
    llpc  = 0.09/0.33;
end
defivalue = llpc*gdp/sum(data.Npop(2:4));
data.vsy  = defivalue;

end
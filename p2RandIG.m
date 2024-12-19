function data = p2RandIG(data,country_data,location)

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
    data.Tres    = lognrnd(3.06, 0.293);
    data.sdl     = 0.1 + 0.9*betarnd(1.71, 2.80);
    data.sdb     = lognrnd(4.03, 2.84);
    data.sdc     = lognrnd(-4.03, 0.873);
    data.t_tit   = wblrnd(75.0, 1.85);
    data.trate   = lognrnd(2.61, 1.47);
    data.Hmax    = lognrnd(2.83, 0.880);
    data.t_vax   = lognrnd(6.32, 0.129);
    data.arate   = wblrnd(167., 1.09);
    data.puptake = wblrnd(0.566, 1.26);
elseif strcmp(location,'UMIC');
    data.Tres    = wblrnd(21.5, 8.40);
    data.sdl     = 0.1 + 0.9*betarnd(1.48, 3.76);
    data.sdb     = lognrnd(4.70, 2.81);
    data.sdc     = lognrnd(-4.16, 0.739);
    data.t_tit   = gamrnd(2.80, 1/0.0702);
    data.trate   = lognrnd(4.47, 0.749);
    data.Hmax    = lognrnd(3.94, 0.793);
    data.t_vax   = gamrnd(94.0, 1/0.199);
    data.arate   = gamrnd(3.53, 1/0.0162);
    data.puptake = wblrnd(0.896, 2.54);
elseif strcmp(location,'HIC');
    data.Tres    = wblrnd(20.2, 8.18);
    data.sdl     = 0.1 + 0.9*betarnd(3.07, 5.32);
    data.sdb     = lognrnd(3.01, 1.66);
    data.sdc     = gamrnd(0.721, 1/89.5);
    data.t_tit   = lognrnd(3.21, 0.584);
    data.trate   = lognrnd(5.62, 0.885);
    data.Hmax    = lognrnd(4.29, 0.738);
    data.t_vax   = lognrnd(6.06, 0.104);
    data.arate   = wblrnd(383, 3.80);
    data.puptake = wblrnd(1.21, 7.08);
elseif any(strcmp(country_data.country,location));
    data.Tres    = country_data{1,startsWith(country_data.Properties.VariableNames, 'Tres')};
    data.sdl     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sdl')};
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
na        = [data.Npop(1:17)',sum(data.Npop(18:end))];%length is 18 to match life table
la        = data.la;
lg        = [dot(la(1),na(1))/sum(na(1)),...
             dot(la(2:4),na(2:4))/sum(na(2:4)),...
             dot(la(5:13),na(5:13))/sum(na(5:13)),...
             dot(la(14:end),na(14:end))/sum(na(14:end))];
lgd       = zeros(size(lg));
for k = 1:length(lgd); 
    lgd(k) = sum(1./((1+0.03).^[1:lg(k)]));
end 
gdp       = 365*sum(data.obj);%in millions USD
gdppc     = gdp/sum(na);
%Masterman & Viscusi (2018) method, using US 2019 VSL of $10.9m from 
%https://www.transportation.gov/office-policy/transportation-policy/revised-departmental-guidance-on-valuation-of-a-statistical-life-in-economic-analysis
if all(strcmp(country_data.igroup,'LLMIC')) || (all(strcmp(country_data.igroup,'UMIC')) && gdppc < 0.008809);
    vsl   = 10.9*((0.008809/0.060362)^0.85)*(gdppc/0.008809);
elseif (all(strcmp(country_data.igroup,'UMIC')) && gdppc > 0.008809) || all(strcmp(country_data.igroup,'HIC'));
    vsl   = 10.9*((gdppc/0.060362)^0.85);
end
defivalue = vsl/(dot(lgd,[na(1);sum(na(2:4));sum(na(5:13));sum(na(14:end))])/sum(na));
data.vly  = defivalue;

%vsy
defivalue = 0.5454*gdp/sum(data.Npop(2:4));
data.vsy  = defivalue;

end
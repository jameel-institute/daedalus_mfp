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
    
    SIG = [...
     1.0000	-0.3461	-0.2443	-0.0367	-0.2517	-0.3313	 0.1709	-0.1753	 0.2787	-0.4842	-0.4689
    -0.3461	 1.0000	 0.2685	 0.0909	 0.1642	 0.4074	 0.1285	 0.2136	-0.2279	 0.4719	 0.4600
    -0.2443	 0.2685	 1.0000	 0.0685	 0.1382	 0.3029	-0.0893	 0.4622	-0.1407	 0.4377	 0.3836
    -0.0367	 0.0909	 0.0685	 1.0000	-0.0592	-0.1391	-0.0463	-0.1191	-0.1213	-0.0379	-0.0967
    -0.2517	 0.1642	 0.1382	-0.0592	 1.0000	 0.4135	 0.1210	 0.2859	-0.3029	 0.1239	 0.2584
    -0.3313	 0.4074	 0.3029	-0.1391	 0.4135	 1.0000	 0.2518	 0.4831	-0.5183	 0.6167	 0.7007
     0.1709	 0.1285	-0.0893	-0.0463	 0.1210	 0.2518	 1.0000	 0.1403	-0.0704	 0.1321	 0.1443
    -0.1753	 0.2136	 0.4622	-0.1191	 0.2859	 0.4831	 0.1403	 1.0000	-0.3580	 0.3235	 0.3137
     0.2787	-0.2279	-0.1407	-0.1213	-0.3029	-0.5183	-0.0704	-0.3580	 1.0000	-0.1678	-0.3863
    -0.4842	 0.4719	 0.4377	-0.0379	 0.1239	 0.6167	 0.1321	 0.3235	-0.1678	 1.0000	 0.8742
    -0.4689	 0.4600	 0.3836	-0.0967	 0.2584	 0.7007	 0.1443	 0.3137	-0.3863	 0.8742	 1.0000];
    
    Z = mvnrnd(zeros(1,11), SIG);
    U = normcdf(Z);

    data.Tres    = logninv(U(1),3.0575,0.2935);
    data.sda     = norminv(U(2),0.2049,1.3057);
    data.sdb     = gaminv(U(3),0.3432,1/0.4712);
    data.sdc     = logninv(U(4),-4.7336,0.8297);
    data.t_tit   = gaminv(U(5),3.1235,1/0.0619);
    data.trate   = logninv(U(6),2.4965,1.5519);
    data.masc    = betainv(U(7),1.2960,1.1635);
    data.Hmax    = logninv(U(8),2.8324,0.8798);
    data.t_vax   = logninv(U(9),6.3189,0.1294);
    data.arate   = wblinv(U(10),168.4494,1.0946);
    data.puptake = betainv(U(11),0.8939,1.6483);

elseif strcmp(location,'UMIC');

    SIG = [...
     1.0000	-0.1835	-0.0783	-0.3563	-0.1713	-0.1970	 0.0324	-0.0640	 0.1534	 0.0131	 0.0716
    -0.1835	 1.0000	 0.3621	 0.4899	-0.1624	-0.1142	 0.0757	-0.0180	-0.1144	 0.3853	 0.2306
    -0.0783	 0.3621	 1.0000	 0.1150	-0.1793	 0.1244	 0.3112	 0.0127	-0.1963	 0.0864	 0.0419
    -0.3563	 0.4899	 0.1150	 1.0000	-0.1431	-0.0790	-0.0723	 0.1449	-0.0208	-0.1115	-0.2484
    -0.1713	-0.1624	-0.1793	-0.1431	 1.0000	-0.0902	 0.1678	-0.0583	-0.0172	-0.3926	-0.3612
    -0.1970	-0.1142	 0.1244	-0.0790	-0.0902	 1.0000	 0.3233	 0.4331	-0.2056	-0.0281	 0.0205
     0.0324	 0.0757	 0.3112	-0.0723	 0.1678	 0.3233	 1.0000	-0.0660	 0.2202	-0.0604	-0.1667
    -0.0640	-0.0180	 0.0127	 0.1449	-0.0583	 0.4331	-0.0660	 1.0000	-0.1968	-0.1998	-0.2294
     0.1534	-0.1144	-0.1963	-0.0208	-0.0172	-0.2056	 0.2202	-0.1968	 1.0000	 0.1207	-0.0592
     0.0131	 0.3853	 0.0864	-0.1115	-0.3926	-0.0281	-0.0604	-0.1998	 0.1207	 1.0000	 0.8732
     0.0716	 0.2306	 0.0419	-0.2484	-0.3612	 0.0205	-0.1667	-0.2294	-0.0592	 0.8732	 1.0000];

    Z = mvnrnd(zeros(1,11), SIG);
    U = normcdf(Z);

    data.Tres    = wblinv(U(1),21.4521,8.3960);
    data.sda     = norminv(U(2),0.6266,1.6400);
    data.sdb     = gaminv(U(3),0.3666,1/0.7220);
    data.sdc     = logninv(U(4),-5.0358,0.6235);
    data.t_tit   = gaminv(U(5),2.8234,1/0.0705);
    data.trate   = logninv(U(6),4.4726,0.7496);
    data.masc    = betainv(U(7),1.5953,1.2964);
    data.Hmax    = logninv(U(8),3.9377,0.7926);
    data.t_vax   = gaminv(U(9),93.9510,1/0.1984);
    data.arate   = gaminv(U(10),3.5593,1/0.0162);
    data.puptake = betainv(U(11),2.2598,2.0971);

elseif strcmp(location,'HIC');
    
    SIG = [...
     1.0000	 0.0750	-0.1962	 0.2533	 0.0856	-0.0156	-0.0941	-0.0969	-0.1407	 0.0053	-0.3129
     0.0750	 1.0000	 0.0999	 0.4183	-0.0714	-0.0788	-0.0535	 0.0701	-0.0240	-0.2640	-0.2102
    -0.1962	 0.0999	 1.0000	-0.1415	-0.4020	 0.0688	-0.4672	 0.4442	-0.0174	 0.1392	-0.0074
     0.2533	 0.4183	-0.1415	 1.0000	-0.3435	 0.1324	 0.1794	 0.0702	-0.3142	-0.1519	-0.1187
     0.0856	-0.0714	-0.4020	-0.3435	 1.0000	-0.2508	 0.0255	 0.0035	 0.2742	-0.0389	-0.1426
    -0.0156	-0.0788	 0.0688	 0.1324	-0.2508	 1.0000	 0.1173	 0.0090	-0.5227	-0.0366	 0.2518
    -0.0941	-0.0535	-0.4672	 0.1794	 0.0255	 0.1173	 1.0000	-0.4514	 0.0311	 0.0653	 0.1012
    -0.0969	 0.0701	 0.4442	 0.0702	 0.0035	 0.0090	-0.4514	 1.0000	-0.0014	-0.2192	-0.3216
    -0.1407	-0.0240	-0.0174	-0.3142	 0.2742	-0.5227	 0.0311	-0.0014	 1.0000	 0.2823	 0.0569
     0.0053	-0.2640	 0.1392	-0.1519	-0.0389	-0.0366	 0.0653	-0.2192	 0.2823	 1.0000	 0.4803
    -0.3129	-0.2102	-0.0074	-0.1187	-0.1426	 0.2518	 0.1012	-0.3216	 0.0569	 0.4803	 1.0000];

    Z = mvnrnd(zeros(1,11), SIG);
    U = normcdf(Z);

    data.Tres    = wblinv(U(1),20.2176,8.1795);
    data.sda     = norminv(U(2),0.2629,0.7784);
    data.sdb     = gaminv(U(3),0.3357,1/0.4750);
    data.sdc     = gaminv(U(4),0.5449,1/141.0678);
    data.t_tit   = logninv(U(5),3.2143,0.5836);
    data.trate   = logninv(U(6),5.6181,0.8854);
    data.masc    = betainv(U(7),1.2295,1.3341);
    data.Hmax    = logninv(U(8),4.2928,0.7379);
    data.t_vax   = logninv(U(9),6.0572,0.1035);
    data.arate   = wblinv(U(10),383.6868,3.7818);
    data.puptake = betainv(U(11),8.1454,2.8712);

elseif any(strcmp(country_data.country,location));

    data.Tres    = country_data{1,startsWith(country_data.Properties.VariableNames, 'Tres')};
    data.sda     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sda')};
    data.sdb     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sdb')};
    data.sdc     = country_data{1,startsWith(country_data.Properties.VariableNames, 'sdc')};
    data.t_tit   = country_data{1,startsWith(country_data.Properties.VariableNames, 't_tit')};
    data.trate   = country_data{1,startsWith(country_data.Properties.VariableNames, 'trate')};
    data.masc    = country_data{1,startsWith(country_data.Properties.VariableNames, 'masc')};
    data.Hmax    = country_data{1,startsWith(country_data.Properties.VariableNames, 'Hmax')};
    data.t_vax   = country_data{1,startsWith(country_data.Properties.VariableNames, 't_vax')};
    data.arate   = country_data{1,startsWith(country_data.Properties.VariableNames, 'arate')};
    data.puptake = country_data{1,startsWith(country_data.Properties.VariableNames, 'puptake')};

    if any(cellfun(@(x) isnumeric(x) && any(isnan(x(:))), struct2cell(data))), error('Not all country data available!'); end

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
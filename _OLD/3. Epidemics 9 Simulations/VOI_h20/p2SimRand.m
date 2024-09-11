function p2SimRand
    
    parpool;
    addpath('../');
    
    igroups    = {'LLMIC','UMIC','HIC'};
    ligroups   = length(igroups);
    load('Argentina.mat','data');%loading Argentina, but only keeping country-independent parameters
    fields     = fieldnames(data);
    ikeep      = [6,7,8,13,14,16,17,18];
    data       = rmfield(data,fields(~ismember(1:numel(fields),ikeep)));  
    lx         = length(data.B);
    data.tvec  = [-75 365*3+1];
    CD         = readtable('country_data.csv');
    nsamples   = 1000;
    diseases   = {'Influenza 2009','Influenza 1918',...
                  'Covid Omicron','Covid Wildtype','SARS'};
    ldis       = length(diseases);
    %strategies = {'No Closures'};%strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
    strategies = {'No Closures','Economic Closures','No Closures','No Closures','Economic Closures';
                  'No Closures','Economic Closures','No Closures','No Closures','Economic Closures';
                  'No Closures','Economic Closures','No Closures','Economic Closures','Economic Closures'};
    lstrat     = 1;%length(strategies);
    inputs     = zeros(nsamples,493);
    output     = zeros(ldis,lstrat,nsamples,2,4);
    %vc         = readtable('40vax_pc_gdp.csv');
    intcost    = zeros(nsamples,1);
    
    for h = 1;%:2;%1:ligroups;
    parfor i = 1:nsamples;
    ldata       = data;
    ldata       = p2RandCountry(ldata,CD,igroups{h});
    [ldata,~,~] = p2Params(ldata,'Covid Wildtype');%to define wnorm and Td_CWT
    inputs(i,:) = [i,NaN,NaN,ldata.Npop',ldata.NNs(1:45)',...
                   ldata.CM(:)',ldata.comm,ldata.travelA3,ldata.schoolA1,ldata.schoolA2,ldata.workp,...
                   ldata.obj',ldata.wfh(1,:),ldata.wfh(2,:),...
                   ldata.t_vax,ldata.arate,ldata.puptake,ldata.Hmax,ldata.t_tit,ldata.trate,ldata.Tres,ldata.sdl,ldata.sdb,...
                   NaN,ldata.la];
    
%     nonempind = find(~isnan(vc.cost_gdp_income40) & strcmp(vc.incomegroup,igroups{h}));
%     randindex = nonempind(randi(numel(nonempind)))   
%     v100g     = (100/40)*vc.cost_gdp_income40(randindex);%cost of vax population as % gdp
%     v100      = v100g*365*sum(ldata.obj); %cost of vax population
%     %units   = 2*sum(ldata.Npop)*(max(ldata.puptake,0.70)-ldata.puptake);%number of vaccines
%     %intcost(i) = units*lognrnd(log(10),log(sqrt(1 + (2/10)^2)))/10^6;%cost per dose for country i in MILLIONS
%     intcost(i)= v100*(max(ldata.puptake,0.70)-ldata.puptake);
    
    hb      = (0.20*ldata.Hmax*sum(ldata.Npop)/10^5);
    hbd     = hb*365;
    gdp     = 365*sum(ldata.obj);
    m       = 799.17;%mean and var of hbd in £
    v       = 535.37^2;
    mu      = log((m^2)/sqrt(v+m^2));
    sigma   = sqrt(log(v/(m^2)+1));
    r       = lognrnd(mu,sigma)*gdp/2176203000000;%cost of hbd in million $
    intcost(i) = hbd*r;
    
    for l = 1:2;  
    if l==2;
        %ldata.puptake = max(ldata.puptake,0.70);
        %ldata.t_vax = min(100,ldata.t_vax);
        ldata.Hmax = 1.2*ldata.Hmax;
    end
               
    for j = 1:ldis;
    inp2           = diseases{j};    
    [ldata,dis,p2] = p2Params(ldata,inp2);
    
    for k = 1:lstrat;
    inp3            = strategies{h,j};
    int             = 5;
    if strcmp(inp3,'Elimination');
        xoptim      = [ones(1*lx,1);ldata.x_econ(:,2);ldata.x_elim(:,1);ones(2*lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(2,:);ldata.wfh(1,:);zeros(2,lx)];
        ldata.imand = [2];
        ldata.inext = [2,2,3,2,5];
    elseif strcmp(inp3,'Economic Closures');
        xoptim      = [ones(2*lx,1);ldata.x_econ(:,2);ldata.x_econ(:,1);ones(lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(inp3,'School Closures');
        xoptim      = [ones(2*lx,1);ldata.x_schc(:,2);ldata.x_schc(:,1);ones(lx,1)];
        ldata.hw    = [zeros(1,lx);ldata.wfh(1,:);ldata.wfh(2,:);ldata.wfh(1,:);zeros(1,lx)];
        ldata.imand = [3];
        ldata.inext = [2,3,3,4,5];
    elseif strcmp(inp3,'No Closures');
        xoptim      = [ones(5*lx,1)];
        ldata.hw    = [zeros(5,lx)];
        ldata.imand = [10];
        ldata.inext = [2,2,5];
    else
        error('Unknown Mitigation Strategy!');
    end
    
    try
        [~,~,g]     = p2Run(ldata,dis,inp3,int,xoptim,p2);
        [cost,~]    = p2Cost(ldata,dis,p2,g);
        sec(1)      = sum(cost([3,6,7:10],:),'all');
        sec(2)      = sum(cost([3],:),'all');
        sec(3)      = sum(cost([6],:),'all');
        sec(4)      = sum(cost([7:10],:),'all');
    catch
        disp([i,j,k]);
        sec         = nan(1,4);
    end
    output(j,k,i,l,:) = sec;
    
    end
    end
    end
    disp(i);
    end
    
    for j = 1:ldis;
    for k = 1:lstrat;
    T                          = array2table([inputs,squeeze(output(j,k,:,1,:)),squeeze(output(j,k,:,2,:)),intcost]);
    T.Properties.VariableNames = {...
    'country','igroup','gnipc',...
    'Npop1','Npop2','Npop3','Npop4','Npop5','Npop6','Npop7','Npop8','Npop9','Npop10',...
    'Npop11','Npop12','Npop13','Npop14','Npop15','Npop16','Npop17','Npop18','Npop19','Npop20','Npop21',...
    'NNs1','NNs2','NNs3','NNs4','NNs5','NNs6','NNs7','NNs8','NNs9','NNs10','NNs11','NNs12','NNs13','NNs14','NNs15',...
    'NNs16','NNs17','NNs18','NNs19','NNs20','NNs21','NNs22','NNs23','NNs24','NNs25','NNs26','NNs27','NNs28','NNs29','NNs30',...
    'NNs31','NNs32','NNs33','NNs34','NNs35','NNs36','NNs37','NNs38','NNs39','NNs40','NNs41','NNs42','NNs43','NNs44','NNs45',...
    'CMaa','CMba','CMca','CMda','CMea','CMfa','CMga','CMha','CMia','CMja','CMka','CMla','CMma','CMna','CMoa','CMpa',...
    'CMab','CMbb','CMcb','CMdb','CMeb','CMfb','CMgb','CMhb','CMib','CMjb','CMkb','CMlb','CMmb','CMnb','CMob','CMpb',...
    'CMac','CMbc','CMcc','CMdc','CMec','CMfc','CMgc','CMhc','CMic','CMjc','CMkc','CMlc','CMmc','CMnc','CMoc','CMpc',...
    'CMad','CMbd','CMcd','CMdd','CMed','CMfd','CMgd','CMhd','CMid','CMjd','CMkd','CMld','CMmd','CMnd','CMod','CMpd',...
    'CMae','CMbe','CMce','CMde','CMee','CMfe','CMge','CMhe','CMie','CMje','CMke','CMle','CMme','CMne','CMoe','CMpe',...
    'CMaf','CMbf','CMcf','CMdf','CMef','CMff','CMgf','CMhf','CMif','CMjf','CMkf','CMlf','CMmf','CMnf','CMof','CMpf',...
    'CMag','CMbg','CMcg','CMdg','CMeg','CMfg','CMgg','CMhg','CMig','CMjg','CMkg','CMlg','CMmg','CMng','CMog','CMpg',...
    'CMah','CMbh','CMch','CMdh','CMeh','CMfh','CMgh','CMhh','CMih','CMjh','CMkh','CMlh','CMmh','CMnh','CMoh','CMph',...
    'CMai','CMbi','CMci','CMdi','CMei','CMfi','CMgi','CMhi','CMii','CMji','CMki','CMli','CMmi','CMni','CMoi','CMpi',...
    'CMaj','CMbj','CMcj','CMdj','CMej','CMfj','CMgj','CMhj','CMij','CMjj','CMkj','CMlj','CMmj','CMnj','CMoj','CMpj',...
    'CMak','CMbk','CMck','CMdk','CMek','CMfk','CMgk','CMhk','CMik','CMjk','CMkk','CMlk','CMmk','CMnk','CMok','CMpk',...
    'CMal','CMbl','CMcl','CMdl','CMel','CMfl','CMgl','CMhl','CMil','CMjl','CMkl','CMll','CMml','CMnl','CMol','CMpl',...
    'CMam','CMbm','CMcm','CMdm','CMem','CMfm','CMgm','CMhm','CMim','CMjm','CMkm','CMlm','CMmm','CMnm','CMom','CMpm',...
    'CMan','CMbn','CMcn','CMdn','CMen','CMfn','CMgn','CMhn','CMin','CMjn','CMkn','CMln','CMmn','CMnn','CMon','CMpn',...
    'CMao','CMbo','CMco','CMdo','CMeo','CMfo','CMgo','CMho','CMio','CMjo','CMko','CMlo','CMmo','CMno','CMoo','CMpo',...
    'CMap','CMbp','CMcp','CMdp','CMep','CMfp','CMgp','CMhp','CMip','CMjp','CMkp','CMlp','CMmp','CMnp','CMop','CMpp',...
    'comm','travelA3','schoolA1','schoolA2','workp',...
    'obj1','obj2','obj3','obj4','obj5','obj6','obj7','obj8','obj9','obj10','obj11','obj12','obj13','obj14','obj15',...
    'obj16','obj17','obj18','obj19','obj20','obj21','obj22','obj23','obj24','obj25','obj26','obj27','obj28','obj29','obj30',...
    'obj31','obj32','obj33','obj34','obj35','obj36','obj37','obj38','obj39','obj40','obj41','obj42','obj43','obj44','obj45',...
    'wfhl1','wfhl2','wfhl3','wfhl4','wfhl5','wfhl6','wfhl7','wfhl8','wfhl9','wfhl10','wfhl11','wfhl12','wfhl13','wfhl14','wfhl15',...
    'wfhl16','wfhl17','wfhl18','wfhl19','wfhl20','wfhl21','wfhl22','wfhl23','wfhl24','wfhl25','wfhl26','wfhl27','wfhl28','wfhl29','wfhl30',...
    'wfhl31','wfhl32','wfhl33','wfhl34','wfhl35','wfhl36','wfhl37','wfhl38','wfhl39','wfhl40','wfhl41','wfhl42','wfhl43','wfhl44','wfhl45',...
    'wfhu1','wfhu2','wfhu3','wfhu4','wfhu5','wfhu6','wfhu7','wfhu8','wfhu9','wfhu10','wfhu11','wfhu12','wfhu13','wfhu14','wfhu15',...
    'wfhu16','wfhu17','wfhu18','wfhu19','wfhu20','wfhu21','wfhu22','wfhu23','wfhu24','wfhu25','wfhu26','wfhu27','wfhu28','wfhu29','wfhu30',...
    'wfhu31','wfhu32','wfhu33','wfhu34','wfhu35','wfhu36','wfhu37','wfhu38','wfhu39','wfhu40','wfhu41','wfhu42','wfhu43','wfhu44','wfhu45',...
    't_vax','arate','puptake','Hmax','t_tit','trate','Tres','sdl','sdb',...
    'vly','la1','la2','la3','la4','la5','la6','la7','la8','la9',...
    'la10','la11','la12','la13','la14','la15','la16','la17','la18',...
    'SEC','VLYL','VSYL','GDPL','SECt','VLYLt','VSYLt','GDPLt','IntCost'};
    writetable(T,strcat('P2L',string(igroups{h}),'_',string(diseases{j}),'_',string(strategies{h,j}),'.csv'));
    end
    end
    end
    %plots = p2Plot(data,f,p2,g,cost,ccost_t,sec(1),inp1,inp2,inp3);
    
    delete(gcp);
    
end
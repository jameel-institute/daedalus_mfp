function [data,Dout,dis,p2]=p2_SARS(data,numInt,inp3)

%data.NNs - column vector of population
%Possible generalsiation to within-sector heterogeneity - one column per subsector

%% COUNTRY PARAMETERS:

%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);

%Population by Sector
ntot = size(data.NNs,1);

%Age-Sector Breakdown
lc=4;
adInd=3;
lx=length(data.NNs)-lc;
[Dout,data]=p2MakeDs(data,data.NNs,ones(lx,1),zeros(1,lx));

%% INITIAL DISEASE PARAMETERS:

dis=struct;%see epidemiological parameter file

%Probabilities
dis.ps = 0.669;
ihr    = dis.ps*8*[0.02284 0.00398 0.00478 0.00983 ...
                 0.01700 0.02922 0.02470 0.02205 ...
                 0.01647 0.01195 0.01647 0.01169 ...
                 0.03081 0.04144 0.04941 0.04941 0.04941];
ifr    = dis.ps*[0.02284 0.00398 0.00478 0.00983 ...
                 0.01700 0.02922 0.02470 0.02205 ...
                 0.01647 0.01195 0.01647 0.01169 ...
                 0.03081 0.04144 0.04941 0.04941 0.04941];
phgs   = ihr./dis.ps;
pdgh   = ifr./ihr;

phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),lx,1);phgs];
nnh     = nn.*ihr;
nnhtot  = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nnhtot  = repelem(nnhtot,ranges);
nnhprop = nnh./nnhtot;
pdgh    = accumarray(subs',pdgh.*nnhprop);
dis.pd  = [repmat(pdgh(adInd),lx,1);pdgh];

%Durations
dis.Tlat  = 1.6;
dis.Tay   = 1.0;
dis.Tsr   = 1.0;
dis.Tsh   = 1.0;
dis.Threc = 5.0;
dis.Thd   = 5.0;
dis.Ti    = 365;

%Calculations
dis.Ts   = ((1-dis.ph).*dis.Tsr)  +(dis.ph.*dis.Tsh);
dis.Th   = ((1-dis.pd).*dis.Threc)+(dis.pd.*dis.Thd);

dis.sig1 = (1-dis.ps)/dis.Tlat;
dis.sig2 = dis.ps/dis.Tlat;
dis.g1   = 1/dis.Tay;%finite for next-generation matrix, pr.p1=1 prevents asymptomatic cases
dis.g2   = (1-dis.ph)./dis.Ts;
dis.g3   = (1-dis.pd)./dis.Th;
dis.h    = dis.ph./dis.Ts;
dis.mu   = dis.pd./dis.Th;
dis.nu   = 1/dis.Ti;

%Transmission
dis.red = 0.58;
dis.R0  = 2.5000;

Deff=Dout.*repmat(data.NNs,1,ntot)./repmat(data.NNs',ntot,1);
onesn=ones(ntot,1);
F=zeros(3*ntot,3*ntot);
F(1:ntot,ntot+1:end)=[dis.red*Deff,Deff];

vvec=[(dis.sig1+dis.sig2).*onesn;      dis.g1.*onesn;       (dis.g2+dis.h).*onesn];%g2 and h are vectors
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-dis.sig1.*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-dis.sig2.*onesn);

GD=F/V;
d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 
dis.beta=dis.R0/R0a;%beta scales actual R0 to fitted R0

%Vaccination
dis.hrv1  = 1/28;                       %time to develop v-acquired immunity
dis.scv1  = 0.60;                       %infection-blocking efficacy
dis.heff  = 0.87;                       %severe-disease-blocking efficacy
dis.hv1   = 1-((1-dis.heff)/(1-dis.scv1)); 
dis.dv1   = 0;                          %death-blocking efficacy
dis.trv1  = 0.52;                       %transmission-blocking efficacy
dis.nuv1  = 1/365;                      %duration of v-acquired immunity

dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
dis.Th_v1 = ((1-(1-dis.dv1)*dis.pd).*dis.Threc)+((1-dis.dv1)*dis.pd.*dis.Thd);

dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
dis.g3_v1 = (1-(1-dis.dv1)*dis.pd)./dis.Th_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;
dis.mu_v1 = (1-dis.dv1)*dis.pd./dis.Th_v1;

%% PREPAREDNESS PARAMETERS:

p2 = struct;%see country parameter file

%Mitigation Time
Hm    = 2*sum(data.Npop)/(10^5);
ihr   = dis.ps*dot(dis.ph,data.NNs)/sum(data.NNs);
Tg    = dis.Tlat+(1-dis.ps)*dis.Tay+dis.ps*((1-(ihr/dis.ps))*dis.Tsr+(ihr/dis.ps)*dis.Tsh);
r     = log(dis.R0)/Tg;
I0    = 10^(-3.22);%
Tm    = (log(Hm/(I0*ihr))/r)-60;
%p2.Tm = Tm-0;
%sector closures (x), working-from-home (wfh), NPIs (delta) and self-isolation (p3/p4) implemented

%Self-Isolation and Quarantine
%p2.p3 = 0.0038;%proportion of asymptomatic self-isolating
%p2.p4 = 0.0038;%proportion of symptomatic self-isolating

%Adherence to NPIs
%pmod       = 0.5801;%lockdown delta%from fitting
%p2.betamod = [1,pmod*ones(1,numInt)];%may be modified in heSingleSim.m or heSwitchSim.m

%Hospital Capacity (per 100k)
%p2.Hmax     = 93.6250*sum(data.Npop)/(10^5);
%p2.thl      = 0.30*p2.Hmax;

dis.Th_oc   = ((1-1.19*dis.pd).*dis.Threc)+(1.19*dis.pd.*dis.Thd);
dis.Th_ocv1 = ((1-1.19*(1-dis.dv1)*dis.pd).*dis.Threc)+(1.19*(1-dis.dv1)*dis.pd.*dis.Thd);

p2.g3_oc    = (1-1.19*dis.pd)./dis.Th_oc;
p2.g3_ocv1  = (1-1.19*(1-dis.dv1)*dis.pd)./dis.Th_ocv1;
p2.mu_oc    = 1.19*dis.pd./dis.Th_oc;%Wilde et al. (2021)
p2.mu_ocv1  = 1.19*(1-dis.dv1)*dis.pd./dis.Th_ocv1;%Wilde et al. (2021)

%Administration Rate (per 100k)
%arate = 682.5457*sum(data.Npop)/(10^5);

%Population
Npop=   data.Npop;
NNage=  [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];

%Uptake
%uptake= 0.7542;

%% PREPAREDNESS LEVELS

if strcmp(inp3,'BAD');
    p2.Tm=Tm+10;
    p2.p3=0;p2.p4=0;
    p2.betamod=[1,0.60*ones(1,numInt)];
    p2.Hmax=12*sum(data.Npop)/(10^5);
    arate=105*sum(data.Npop)/(10^5);
    uptake= 0.60;
    
elseif strcmp(inp3,'MEDIUM');
    p2.Tm=Tm-10;
    p2.p3=0.0125;p2.p4=0.0125;
    p2.betamod=[1,0.50*ones(1,numInt)];
    p2.Hmax=132*sum(data.Npop)/(10^5);
    arate=657.5*sum(data.Npop)/(10^5);
    uptake= 0.775;

elseif strcmp(inp3,'GOOD');    
    p2.Tm=Tm-30;
    p2.p3=0.025;p2.p4=0.025;
    p2.betamod=[1,0.40*ones(1,numInt)];
    p2.Hmax=252*sum(data.Npop)/(10^5);
    arate=1210*sum(data.Npop)/(10^5);
    uptake= 0.95;
end

p2.thl=0.30*p2.Hmax;

%%

%Rollout
t_vax=  366;%1st January (of second year)
%tpoints=cumsum([t_vax,uptake*(fliplr(NNage)./arate)]);
t_ages= uptake*(NNage./arate);
t_ages= [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
tpoints=cumsum([t_vax,t_ages]);

%Period 1 - retired-age
p2.startp1= tpoints(1);
p2.aratep1= [0;0;0;arate];

%Period 2 - working-age
p2.startp2= tpoints(2);
p2.aratep2= [0;0;arate;0];%to be split across all economic sectors in heSimCovid19vax.m

%Period 3 - school-age
p2.startp3= tpoints(3);
p2.aratep3= [0;arate;0;0];

%Period 4 - pre-school-age
p2.startp4= tpoints(4);
p2.aratep4= [arate;0;0;0];

%End of Rollout
p2.end=     tpoints(5);

end
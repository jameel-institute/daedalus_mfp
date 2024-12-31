function [dis,p2] = dd_set_disease(data,inp2)

%% COUNTRY PARAMETERS:

%Population by Age
nn     = data.Npop';
nn     = [nn(1:16),sum(nn(17:end))];%last age range for disease is 80+
nntot  = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
ranges = [1,3,9,4];
nntot  = repelem(nntot,ranges);
nnprop = nn./nntot;
subs   = 1:4;
subs   = repelem(subs,ranges);

%Population by Sector
ln          = length(data.NNs);
lx          = length(data.obj);
adInd       = 3;

%Contact Matrix
Dout = dd_calc_contacts(data,ones(lx,1),zeros(1,lx));

%% INITIAL DISEASE PARAMETERS:

addpath('input');

if strcmp(inp2,'Influenza 2009');
    dis = param_influenza_2009;
elseif strcmp(inp2,'Influenza 1957');
    dis = param_influenza_1957;
elseif strcmp(inp2,'Influenza 1918');
    dis = param_influenza_1918;
elseif strcmp(inp2,'Covid Omicron');
    dis = param_covid_omicron;    
elseif strcmp(inp2,'Covid Delta');
    dis = param_covid_delta;   
elseif strcmp(inp2,'Covid Wildtype');
    dis = param_covid_wildtype;    
elseif strcmp(inp2,'SARS');
    dis = param_sars;
else
    error('Unknown Disease!');
end   

%Probabilities
phgs    = dis.ihr./dis.ps;
pdgh    = dis.ifr./dis.ihr;
phgs    = accumarray(subs',phgs.*nnprop);
dis.ph  = [repmat(phgs(adInd),lx,1);phgs];
nnh     = nn.*dis.ihr;
nnhtot  = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nnhtot  = repelem(nnhtot,ranges);
nnhprop = nnh./nnhtot;
pdgh    = accumarray(subs',pdgh.*nnhprop);
dis.pd  = [repmat(pdgh(adInd),lx,1);pdgh];

%Durations
dis.Ts = ((1-dis.ph).*dis.Tsr)   + (dis.ph.*dis.Tsh);
dis.Th = ((1-dis.pd).*dis.Threc) + (dis.pd.*dis.Thd);

%Rates
dis.siga = (1-dis.ps)/dis.Tlat;
dis.sigs = dis.ps/dis.Tlat;
dis.g1   = 1/dis.Tay;
dis.g2   = (1-dis.ph)./dis.Ts;
dis.g3   = (1-dis.pd)./dis.Th;
dis.h    = dis.ph./dis.Ts;
dis.mu   = dis.pd./dis.Th;
dis.nu   = 1/dis.Ti;

%Vaccination
dis.hrv1 = 1/21;%1/28;                 %time to develop v-acquired immunity
dis.scv1 = 0.70;%0.60;                 %infection-blocking efficacy
dis.heff = 0.85;%0.87;                 %severe-disease-blocking efficacy
dis.hv1  = 1-((1-dis.heff)/(1-dis.scv1)); 
dis.trv1 = 0.30;%0.52;                 %transmission-blocking efficacy
dis.nuv1 = 1/365;                      %duration of v-acquired immunity

dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr)  +((1-dis.hv1)*dis.ph.*dis.Tsh);
dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;

%Transmission
%dis.beta = 1;
[r0,ev]   = dd_calc_Rt(dis,dis.h,dis.g2,data.NNs,zeros(ln,1),zeros(ln,1),data.NNs,Dout,1,dis.siga,dis.sigs,0,0,1,1);
%dis.r0a  = r0;
%dis.beta = dis.R0/dis.r0a;
%dis.r0   = dis.R0;
dis.r0    = r0;
ev        = abs(ev(1:ln));%corresponding eigenvector to seed exposed population
dis.Ev    = ev./sum(ev);

%% PREPAREDNESS PARAMETERS:

p2 = struct;

p2.Tres  = data.Tres;                       %Response Time
p2.t_tit = data.t_tit;                      %Test-Isolate-Trace Time
p2.trate = data.trate;                      %Test-Isolate-Trace Rate
p2.sda   = data.sda;                        %Distancing Multiplier Minimum
p2.sdb   = data.sdb;                        %Distancing Death Sensitivity
p2.sdc   = data.sdc;                        %Distancing Time Relaxation
p2.Hmax  = data.Hmax*sum(data.Npop)/10^5;   %Hospital Capacity
t_vax    = data.t_vax;                      %Vaccine Administration Time
arate    = data.arate*sum(data.Npop/10^5);  %Vaccine Administration Rate
puptake  = data.puptake;                    %Vaccine Uptake

%Response Time
r       = dd_calc_r(dis,dis.h,dis.g2,dis.mu,dis.g3,data.NNs, ...
                    zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1), ...
                    data.NNs,Dout,1,dis.siga,dis.sigs,0,0,1,1);
Td      = log(2)/r;
p2.Tres = p2.Tres*Td;

%Test-Isolate-Trace
p2.t_tit = p2.t_tit*Td;

%Hospital Capacity

%Vaccine Uptake
Npop    = data.Npop;
NNage   = [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];
puptake = puptake*(1-(1/dis.r0));
puptake = min(puptake,0.95*(1-(NNage(1)/sum(NNage))));%population uptake cannot be greater than 95% coverage in non-pre-school age groups
upfun   = @(up) puptake*sum(NNage) - min(0.5*up,0.95)*NNage(2) - min(up,0.95)*NNage(3) - min(1.5*up,0.95)*NNage(4);
try
    up  = fzero(upfun,[0 2]);
catch
    up  = fminbnd(upfun,0,2);
end
u1      = 0;
u2      = min(0.5*up,0.95);
u3      = min(up,0.95);
u4      = min(1.5*up,0.95);
uptake  = [u1,u2,u3,u4];

%Vaccine Administration Rate
t_ages     = min((uptake.*NNage)/arate,Inf);%arate may be 0
if strcmp(inp2,'Influenza 1918');
    t_ages  = [t_ages(3),t_ages(4),t_ages(2),t_ages(1)];
    aratep1 = [0;0;arate;0];%Period 1 - working-age%to be split across all economic sectors in heSimCovid19vax.m
    aratep2 = [0;0;0;arate];%Period 2 - retired-age
    aratep3 = [0;arate;0;0];%Period 3 - school-age
    aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
else
    t_ages  = [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
    aratep1 = [0;0;0;arate];%Period 1 - retired-age
    aratep2 = [0;0;arate;0];%Period 2 - working-age%to be split across all economic sectors in heSimCovid19vax.m
    aratep3 = [0;arate;0;0];%Period 3 - school-age
    aratep4 = [0;0;0;0];    %Period 4 - pre-school-age
end
NNprop                  = ones(ln,1);
workingage              = sum(data.NNs([1:lx,lx+adInd]));
NNprop([1:lx,lx+adInd]) = data.NNs([1:lx,lx+adInd])/workingage;
p2.aratep1              = NNprop.*[repmat(aratep1(3),lx,1);aratep1];    
p2.aratep2              = NNprop.*[repmat(aratep2(3),lx,1);aratep2];
p2.aratep3              = NNprop.*[repmat(aratep3(3),lx,1);aratep3];
p2.aratep4              = NNprop.*[repmat(aratep4(3),lx,1);aratep4];

%Vaccine Administration Time
tpoints    = cumsum([t_vax,t_ages]);
p2.startp1 = tpoints(1);
p2.startp2 = tpoints(2);
p2.startp3 = tpoints(3);
p2.startp4 = tpoints(4);
p2.end     = tpoints(5);%End of Rollout

%% COST PARAMETERS:

na     = [data.Npop(1:16)',sum(data.Npop(17:end))];%length is 17 to match ifr
la     = [data.la(1:16),...
          dot(data.la(17:end),[data.Npop(17),sum(data.Npop(18:end))])/sum(data.Npop(17:end))];
napd   = na.*dis.ifr;
lg     = [dot(la(1),napd(1))/sum(napd(1)),...
          dot(la(2:4),napd(2:4))/sum(napd(2:4)),...
          dot(la(5:13),napd(5:13))/sum(napd(5:13)),...
          dot(la(14:end),napd(14:end))/sum(napd(14:end))];
dis.lg = [repmat(lg(3),1,45),lg];

end
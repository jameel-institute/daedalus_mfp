function [dis,p2] = dd_set_disease(data,inp2)

ln    = length(data.NNs);
lx    = length(data.obj);
adInd = 3;

%% DISEASE PARAMETERS:

%addpath('input');
dis = feval(strcat('param_',lower(strrep(inp2,' ','_'))));

%Population by Age
no      = data.Npop';
nn      = [no(1:16),sum(no(17:end))];%last age range for disease is 80+
nnh     = nn.*dis.ihr;
nnd     = nn.*dis.ifr;
nnc     = [nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
nnhc    = [nnh(1),sum(nnh(2:4)),sum(nnh(5:13)),sum(nnh(14:end))];
nndc    = [nnd(1),sum(nnd(2:4)),sum(nnd(5:13)),sum(nnd(14:end))];
nntot   = repelem(nnc,[1,3,9,4]);
nnhtot  = repelem(nnhc,[1,3,9,4]);
nndtot  = repelem(nndc,[1,3,9,4]);
nnprop  = nn./nntot;
nnhprop = nnh./nnhtot;
nndprop = nnd./nndtot;
subs    = [1:4];
subs    = repelem(subs,[1,3,9,4]);

%Proportions
phgs   = dis.ihr./dis.ps;
pdgh   = dis.ifr./dis.ihr;
phgs   = accumarray(subs',phgs.*nnprop);
pdgh   = accumarray(subs',pdgh.*nnhprop);
dis.ph = [repmat(phgs(adInd),lx,1);phgs];
dis.pd = [repmat(pdgh(adInd),lx,1);pdgh];

%Durations
dis.Ts = ((1-dis.ph).*dis.Tsr)   + (dis.ph.*dis.Tsh);
dis.Th = ((1-dis.pd).*dis.Threc) + (dis.pd.*dis.Thd);

%Transition Rates
dis.siga = (1-dis.ps)/dis.Tlat;
dis.sigs = dis.ps/dis.Tlat;
dis.g1   = 1/dis.Tay;
dis.g2   = (1-dis.ph)./dis.Ts;
dis.g3   = (1-dis.pd)./dis.Th;
dis.h    = dis.ph./dis.Ts;
dis.mu   = dis.pd./dis.Th;
dis.nu   = 1/dis.Ti;

%Vaccine Proportions
dis.scv1 = 0.70;%infection-blocking efficacy
dis.heff = 0.90;%severe-disease-blocking efficacy
dis.trv1 = 0.30;%transmission-blocking efficacy
dis.hv1  = 1-((1-dis.heff)/(1-dis.scv1)); 

%Vaccine Durations
Tsc       = 21; %time to develop v-acquired immunity
dis.Ts_v1 = ((1-(1-dis.hv1)*dis.ph).*dis.Tsr) + ((1-dis.hv1)*dis.ph.*dis.Tsh);
Tiv       = 365;%duration of v-acquired immunity

%Vaccine Transition Rates
dis.hrv1  = 1/Tsc;                       
dis.g2_v1 = (1-(1-dis.hv1)*dis.ph)./dis.Ts_v1;
dis.h_v1  = (1-dis.hv1)*dis.ph./dis.Ts_v1;
dis.nuv1  = 1/Tiv;                      

%Transmission
Dout      = dd_calc_contacts(data,ones(lx,1),zeros(1,lx));
%dis.beta = 1;
[r0,ev]   = dd_calc_Rt(dis,dis.h,dis.g2,data.NNs,zeros(ln,1),zeros(ln,1),data.NNs,Dout,1,dis.siga,dis.sigs,0,0,1,1);
%dis.r0a  = r0;
%dis.beta = dis.R0/dis.r0a;
%dis.r0   = dis.R0;
dis.r0    = r0;
ev        = abs(ev(1:ln));%corresponding eigenvector to seed exposed population
dis.Ev    = ev./sum(ev);

%Initial Condition
zn     = zeros(ln,1);
R0     = zn;%entirely naive population
y0     = [data.NNs-R0;repmat(zn,6,1);R0;repmat(zn,11,1);data.NNs-R0];
dis.y0 = y0;%seeding is defined in dd_run_sim.m

%% OTHER DISEASE-DEPENDENT PARAMETERS:

if nn(end) ~= 0;
    la = [data.la(1:16), dot(data.la(17:end),[no(17),sum(no(18:end))])/nn(end)];
else
    la = [data.la(1:16), 0];
end
lg     = accumarray(subs',la.*nndprop)';
dis.lg = [repmat(lg(adInd),1,lx),lg];

p2 = struct;

%Distancing
Tres    = data.Tres;%response time (in terms of doubling times)
p2.sda  = data.sda; %distancing multiplier intercept
p2.sdb  = data.sdb; %distancing death sensitivity
p2.sdc  = data.sdc; %distancing time relaxation
r       = dd_calc_r(dis,dis.h,dis.g2,dis.mu,dis.g3,data.NNs, ...
                    zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1),zeros(ln,1), ...
                    data.NNs,Dout,1,dis.siga,dis.sigs,0,0,1,1);
Td      = log(2)/r;
p2.Tres = Tres*Td;

%Surveillance
t_tit    = data.t_tit;%case-isolation-tracing start time (in terms of doubling times)
p2.trate = data.trate;%case-isolation-tracing testing rate
p2.t_tit = t_tit*Td;
p2.asca  = 0.7623;%ascertainment proportion coefficients
p2.ascb  = 1.605;
p2.ascc  = -1.416;
p2.pcta  = 2.159;%ascertainment method coefficients
p2.pctb  = 1.697;
p2.opsa  = 11.3224;%symptom-driven onset-PCR delay coefficients  
p2.opsb  = -2.6260;
p2.opc   = -5.6304;%traced onset-PCR delay coefficient

%Healthcare
p2.Hmax = data.Hmax*sum(nn/10^5);%spare hospital capacity
p2.th   = 1.87;%Johnson et al. (2023)

%Vaccination
t_vax   = data.t_vax;             %vaccine administration start time
arate   = data.arate*sum(nn/10^5);%vaccine administration rate
puptake = data.puptake;           %vaccine uptake
%Uptake
av_ifr  = dot(dis.ifr,nn)/sum(nn);
puptake = puptake*(1-(1/r0))*(4^log10(100*av_ifr));
puptake = min(puptake,0.95*(1-(nnc(1)/sum(nnc))));%population uptake cannot be greater than 95% coverage in non-pre-school age groups
upfun   = @(up) puptake*sum(nnc) - min(0.5*up,0.95)*nnc(2) - min(up,0.95)*nnc(3) - min(1.5*up,0.95)*nnc(4);
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
%Administration Rate
t_ages      = min((uptake.*nnc)/arate,Inf);%prevent NaN if numerator and denominator both zero
if strcmp(inp2,'Influenza 1918');
    t_ages  = [t_ages(3),t_ages(4),t_ages(2),t_ages(1)];
    aratep1 = [0;0;arate;0];%period 1 - working-age
    aratep2 = [0;0;0;arate];%period 2 - retired-age
    aratep3 = [0;arate;0;0];%period 3 - school-age
    aratep4 = [0;0;0;0];    %period 4 - pre-school-age
else
    t_ages  = [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
    aratep1 = [0;0;0;arate];%period 1 - retired-age
    aratep2 = [0;0;arate;0];%period 2 - working-age
    aratep3 = [0;arate;0;0];%period 3 - school-age
    aratep4 = [0;0;0;0];    %period 4 - pre-school-age
end
nadprop                  = ones(ln,1);
nadprop([1:lx,lx+adInd]) = data.NNs([1:lx,lx+adInd])/sum(data.NNs([1:lx,lx+adInd]));%adult population proportion
p2.aratep1               = nadprop.*[repmat(aratep1(adInd),lx,1);aratep1];    
p2.aratep2               = nadprop.*[repmat(aratep2(adInd),lx,1);aratep2];
p2.aratep3               = nadprop.*[repmat(aratep3(adInd),lx,1);aratep3];
p2.aratep4               = nadprop.*[repmat(aratep4(adInd),lx,1);aratep4];
%Administration Time
tpoints    = cumsum([t_vax,t_ages]);
p2.startp1 = tpoints(1);
p2.startp2 = tpoints(2);
p2.startp3 = tpoints(3);
p2.startp4 = tpoints(4);
p2.end     = tpoints(5);%end of rollout

end
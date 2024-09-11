function [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepSwineFlu(data,numInt,inp1,inp2)%,R0,del)

%data.NNs - column vector of population
%Possible generalsiation to within-sector heterogeneity - one column per subsector

%% POPULATION PARAMETERS:

%Population Density
[n,na]=size(data.NNs);
ntot=n*na;
NN=sum(data.NNs,2);
NNbar=reshape(data.NNs,ntot,1);
NNrep=repmat(NN,na,1);

%Urban and Rural?
%urbrur=0;%Turn in to @home vs @work7/291
%{
if urbrur==1
    hvec=kron(hvec,ones(n,1));
    muvec=kron(muvec,ones(n,1));
    %Kkron=[1,.05;.05,.75];%Sonoma
    Kkron=[1,.05;.05,.85];%Mendocino
else
    Kkron=1;
end
%}
%{
C=[.3827    0.9115    0.0419;
    1.2062    7.3538    0.1235;
    0.1459    0.4810    0.1860];
%}
%C=eye(na);

%Age-Sector Breakdown
lc=4;
adInd=3;
lx=length(data.NNs)-lc;
D=heMakeDs(NN,ones(lx,1),data,zeros(1,lx));
Dout=D;

%K=heMakeDs(NN,eye(10));
%K=rand(n);
%K=normr(K);
%D=kron(C,K);

%% DISEASE PARAMETERS:

pr=struct;%see epidemiological parameter file

%Transmission
pr.R0=1.7000;
pr.red=0.50;

%Latency and Onset
Text=1.6;
pr.sigma=1/Text;
%Tonset=1;
%pr.omega=1/Tonset;

%Case Pathways
pr.p1=0.669;
[ph,pd,~,~]=heParamsAge(data,pr.p1);
% pili=[0.4670,0.4786,0.6590,0.7252]';
% pili=0.505*pili;
% pili=[repmat(pili(adInd),lx,1);pili];
% pr.p2=pili;
% ph=[0.0500    0.0022    0.0350    0.5235]';
ph=[repmat(ph(adInd),lx,1);ph];
% pd=[0.0103    0.0078    0.0361    0.1555]';
pd=[repmat(pd(adInd),lx,1);pd];
% pdeath=0.39;
% pd=pd/sum(pd);
% pd=pdeath*pd;%new data: 39% of hospital cases result in death

%Hospitalisation and Death Rates (proportion and rate combined)
Tsh=2;
pr.h=ph/Tsh;
Thd=5;
pr.mu=pd/Thd;

%Recovery Rates (proportion and rate combined)
Ta=2;%asymptomatic
pr.g1=1/Ta;
% pr.gX=1/2.1;%mild
Ts=2;%symptomatic
pr.g2=(1-ph)/Ts;
Threc=5;%hospitalised
pr.g3=(1-pd)/Threc;

%Immunity Loss
Ti=inf;
pr.nu=1/Ti;

%%

Deff=Dout.*repmat(NNbar,1,n*na)./repmat(NNbar',n*na,1);
onesn=ones(ntot,1);

F=zeros(3*ntot,3*ntot);
%F=zeros(4*ntot,4*ntot);
F(1:ntot,ntot+1:end)=[pr.red*Deff,Deff];
%F(1:ntot,ntot+1:end)=[pr.red*Deff,repmat(Deff,1,2)];

vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       (pr.g2+pr.h).*onesn];%g2 and h are vectors
%vvec=[pr.sigma.*onesn;      pr.g1.*onesn;       pr.g2.*onesn;       (pr.gX+pr.h).*onesn];%gX and h are vectors
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=    diag(-(1-pr.p1) .*pr.sigma  .*onesn);
V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*pr.sigma  .*onesn);
%V(2*ntot+1:3*ntot,1:ntot)=  diag(-pr.p1     .*(1-pr.p2) .*pr.sigma  .*onesn);
%V(3*ntot+1:4*ntot,1:ntot)=  diag(-pr.p1     .*pr.p2     .*pr.sigma  .*onesn);

GD=F/V;

%Ceff=kron(C,Ckron);%Urb/rural mixing here
%F(1:ntot,1:ntot)=Deff;
%vvec=kron([pr.sigma;pr.g1;pr.g2;pr.gX],ones(ntot,1));
%vvec(end-ntot+1:end)=vvec(end-ntot+1:end)+pr.h;

%{
%HE:
ntot=n*na;
F=zeros(7*ntot,7*ntot);
F(1:ntot,1:ntot)=Deff;%ntot+1:end)=[repmat(2/3*Deff,1,2),repmat(Deff,1,4)];
onesn=ones(ntot,1);
vvec=[pr.sigma*onesn;pr.g1*onesn;pr.omega*onesn;pr.g2*onesn;pr.g2*onesn;pr.h*onesn;pr.h*onesn];
V=diag(vvec);
V(ntot+1:2*ntot,1:ntot)=diag(-pr.sigma*(1-pr.p1)*onesn);
V(2*ntot+1:3*ntot,1:ntot)=diag(-pr.sigma*pr.p1*onesn);
V(3*ntot+1:4*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p3).*(1-pr.p2));
V(4*ntot+1:5*ntot,2*ntot+1:3*ntot)=diag(-pr.p3.*(1-pr.p2));
V(5*ntot+1:6*ntot,2*ntot+1:3*ntot)=diag(-(1-pr.p4).*pr.p2);
V(6*ntot+1:7*ntot,2*ntot+1:3*ntot)=diag(-pr.p4.*pr.p2);
GD=F/V;
%}

d=eigs(GD,1);%largest in magnitude (+/-) 
R0a=max(d); 
beta=pr.R0/R0a;%beta scales actual R0 to fitted R0

%% PREPAREDNESS PARAMETERS:

%see country parameter file

%Mitigation Time
Hm=2*sum(data.Npop)/(10^5);
ihr=pr.p1*dot(ph,NN)/sum(NN);
Tg=Text+(1-pr.p1)*Ta+pr.p1*((1-(ihr/pr.p1))*Ts+(ihr/pr.p1)*Tsh);
r=log(pr.R0)/Tg;
I0=10^(-3.60);%
Tm=(log(Hm/(I0*ihr))/r)-60;
pr.Tm=Tm-0;
%sector closures (x), working-from-home (wfh), NPIs (delta) and self-isolation (p3/p4) implemented

%Adherence to NPIs
pmod=0.5267;%lockdown delta%from fitting
% pmodIn=1.30*pmod;%post-lockdown delta
pr.betamod=[1,pmod*ones(1,numInt)];%may be modified in heSingleSim.m or heSwitchSim.m

%Self-Isolation and Quarantine
pr.p3=0.0003;%proportion of asymptomatic self-isolating
pr.p4=0.0003;%proportion of symptomatic self-isolating
%pr.p5=0.00;%proportion of severe self-isolating
%
pr.odds=0;%asymptomatic quarantining rate
pr.q1=0;%mild quarantining rate
pr.q2=0;%severe quarantining rate
% pr.g4=1/(1/pr.g2-1/pr.q1);%asymptomatic/mild quarantining recovery rate
% pr.qnew=0;
% if pr.q2>0
%     pr.g4X=1/(1./pr.gX+1./pr.q2);%Vector
% else
%     pr.g4X=pr.q2*0;
% end

%Hospital Capacity (per 100k)
pr.Hmax=66.4339*sum(data.Npop)/(10^5);
%pr.thu=0.90*pr.Hmax;
pr.thl=0.30*pr.Hmax;
pr.mu_oc=1.19*pr.mu;%Wilde et al. (2021)
pr.g3_oc=(1-1.19*pd)/Threc;
%see over-capacity death and recovery rates for vaccinated individuals below

%% VACCINATION PARAMETERS:

vx=struct;

%Vaccine 1
vx.hrv1=    1/21;                       %time to develop v-acquired immunity (AstraZeneca)
vx.scv1=    0.61;                       %infection-blocking efficacy
vx.p1v1=    0;                          %disease-blocking efficacy          
vx.hv1=     0;%1-((1-0.61)/(1-vx.scv1));%severe-disease-blocking efficacy
vx.dv1=     0;                          %death-blocking efficacy
vx.trv1=    0;                          %transmission-blocking efficacy
vx.nuv1=    1/Inf;                      %duration of v-acquired immunity

vx.h_v1=    (1-vx.hv1)*ph/Tsh;
vx.g2_v1=   (1-(1-vx.hv1)*ph)/Ts;
vx.mu_v1=   (1-vx.dv1)*pd/Thd;
vx.g3_v1=   (1-(1-vx.dv1)*pd)/Threc;

vx.mu_ocv1= 1.19*vx.mu_v1;%Wilde et al. (2021)
vx.g3_ocv1= (1-1.19*(1-vx.dv1)*pd)/Threc;

% %Vaccine 2
% vx.hrv2=    1/21;                       %time to develop v-acquired immunity (AstraZeneca)
% vx.scv2=    0.65;                       %infection-blocking efficacy
% vx.p1v2=    0;                          %disease-blocking efficacy          
% vx.hv2=     1-((1-0.80)/(1-vx.scv1));   %severe-disease-blocking efficacy
% vx.dv2=     0;                          %death-blocking efficacy
% vx.trv2=    0;                          %transmission-blocking efficacy
% vx.nuv2=    1/Inf;                      %duration of v-acquired immunity

%Administration Rate (per 100k)
arate=152.3081*sum(data.Npop)/(10^5);

%Population
Npop=   data.Npop;
NNage=  [Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:end))];

%Uptake
uptake= 0.8862;

%Rollout
t_vax=  366;%1st January (of second year)
%tpoints=cumsum([t_vax,uptake*(fliplr(NNage)./arate)]);
t_ages= uptake*(NNage./arate);
t_ages= [t_ages(4),t_ages(3),t_ages(2),t_ages(1)];
tpoints=cumsum([t_vax,t_ages]);

%Period 1 - retired-age
vx.startp1= tpoints(1);
vx.aratep1= [0;0;0;arate];%to be split across all economic sectors in heSimCovid19vax.m

%Period 2 - working-age
vx.startp2= tpoints(2);
vx.aratep2= [0;0;arate;0];

%Period 3 - school-age
vx.startp3= tpoints(3);
vx.aratep3= [0;arate;0;0];

%Period 4 - pre-school-age
vx.startp4= tpoints(4);
vx.aratep4= [arate;0;0;0];

%End of Rollout
vx.end=     tpoints(5);

end

%%

function [phgs,pdgh,Threc,Thd]=heParamsAge(datax,ps)
%%

%nn=[3463,3726,3538,3260,3693,4022,4011,3926,3586,3919,4129,3890,3308,2982,2960,2069,1531+933+414+117+12];%England and Wales, mid-2019 estimate (ONS)
nn=datax.Npop';
nn=[nn(1:16),sum(nn(17:end))];%last age range in Knock et al. (2020) is 80+

ranges=[1,3,9,4];
nntot=[nn(1),sum(nn(2:4)),sum(nn(5:13)),sum(nn(14:end))];
nntot=repelem(nntot,ranges);
nnprop=nn./nntot;

subs=1:4;
subs=repelem(subs,ranges);

%%

% ps=     0.6;
ihr=    ps*[0.00697	0.00274	0.00274	0.00274	0.00274	0.00561	0.00561	0.00561	0.00561	0.00561	0.01060	0.01060	0.01060	0.01546	0.01546	0.01546	0.01546];
phgs=   ihr./ps;

picu=   [NaN];
pg=     1-picu;
pdicu=  [NaN];
psd=    1-pdicu;

ifr=    ps*[0.00028	0.00011	0.00011	0.00012	0.00012	0.00030	0.00030	0.00030	0.00030	0.00065	0.00065	0.00065	0.00065	0.00980	0.00980	0.00980	0.00980];
pdgh=   ifr./ihr;

phgs=   accumarray(subs',phgs.*nnprop);
pdgh=   accumarray(subs',pdgh.*nnprop);

%%

Tgr=    10.7;

Ttr=    2.5; 
Ticur=  15.6; 
Tsdr=   12.2;

Tgd=    10.3;

%Ttr=    2.5; 
Ticud=  11.8;

%Ttr=    2.5; 
Ticusd= 7; 
Tsdd=   8.1;

Threc=  Tgr*(pg*nn'/sum(nn))+(Ttr+Ticur+Tsdr)*(picu*nn'/sum(nn));
Thd=    Tgd*(pg*nn'/sum(nn))+(Ttr+Ticud)*((picu.*pdicu)*nn'/sum(nn))+(Ttr+Ticusd+Tsdd)*((picu.*psd)*nn'/sum(nn));

end
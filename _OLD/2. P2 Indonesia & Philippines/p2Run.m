function [tvec,f,g]=p2Run(data,Dout,tvec,dis,Xit,p2)

lt=length(tvec);

numseed=7;

lc    = 4;%*age
adInd = 3;%Age group of adults *age
EdInd = 41;%education sector index
lx    = length(data.obj);%Number of sectors

NNbar=data.NNs;
XitMat=reshape(Xit,lx,lt-2);
WitMat=XitMat.^(1/data.alp);
WitMat(EdInd,:)=XitMat(EdInd,:);
NNvec=repmat(NNbar(1:lx),1,lt-2).*WitMat;%Assumes pre-lockdown=fully open
NNworkSum=sum(NNvec,1);
NNvec(lx+1:lx+lc,1:lt-2)=repmat(NNbar(lx+1:lx+lc),1,lt-2);
NNvec(lx+adInd,:)=sum(NNbar([1:lx,lx+adInd]))-NNworkSum;
NNvec=[NNbar,NNvec];

Dvec=repmat(Dout,[1,1,lt-1]);

for i=2:lt-1
    [Dtemp,~]   = p2MakeDs(data,NNvec(:,i),XitMat(:,i-1),data.wfh(i,:));
    Dvec(:,:,i) = p2.betamod(i)*Dtemp;
end

seed=10^-(numseed);%*NNprob;
seedvec=seed*ones(lx+lc,1); %zeros(ntot,1); seedvec(2*n+1:3*n)=seed*ones(n,1);

[tvec,f,g]=p2SimVax(data,NNvec,Dvec,tvec,dis,NNvec(:,1),seedvec,p2);

end

%%

function [tvec,f,g]=p2SimVax(data,NNvec,Dvec,tvec,dis,S0,seedvec,p2)               
%% PARAMETERS:
ntot=size(data.NNs,1);
adInd=3;
lx=ntot-4;
NNbar=NNvec(:,1);
sumWorkingAge=sum(NNbar([1:lx,lx+3]));

nc=20;
hospInc=0;
    
lt=length(tvec);
zn=zeros(ntot,1);

t0=tvec(1);
y0=[S0;repmat(zn,6,1);NNbar-S0;repmat(zn,nc-9,1);S0];

toutAll=t0;
Sout=S0';
Sv1out=repmat(zn',1,2);
Iout=sum(seedvec');
Hout=zn';
HnewAll=[];
Dout=zn';
DEout=repmat(zn,1,lt);
%Rout=sum((NNbar-S0)');
Vout=zn';
Snout=S0';
Rt=zeros(lt-1,1);

%% FIXED INPUT:

for i=1:lt-1

    tend=tvec(i+1);
    NNfeed=NNvec(:,i);
    NNfeed(NNfeed==0)=1;
    D=Dvec(:,:,i);

    %Vaccination Rollout by Sector
    NNnext=NNvec(:,i);
    NNnext(lx+[1,2])=1;
    NNnext([1:lx,lx+3])=NNnext([1:lx,lx+3])/sumWorkingAge;
    NNnext(end)=1;
    p2.ratep1=NNnext.*[repmat(p2.aratep1(3),lx,1);p2.aratep1];    
    p2.ratep2=NNnext.*[repmat(p2.aratep2(3),lx,1);p2.aratep2];
    p2.ratep3=NNnext.*[repmat(p2.aratep3(3),lx,1);p2.aratep3];
    p2.ratep4=NNnext.*[repmat(p2.aratep4(3),lx,1);p2.aratep4];
    %p2.ratep5=NNnext.*[repmat(p2.aratep5(3),lx,1);p2.aratep5];

    [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0,Sv1class,Vclass,Snclass]=integr8(data,NNfeed,D,i,t0,tend,dis,seedvec,y0,p2);

    toutAll=[toutAll;tout(2:end)];
    Sout=[Sout;Sclass(2:end,:)];
    Sv1out=[Sv1out;Sv1class(2:end,:)];
    Iout=[Iout;Itot(2:end)];
    Hout=[Hout;Hclass(2:end,:)];
    Dout=[Dout;Dclass(2:end,:)];
    DEout(:,i)=DEcum;
    %Rout(:,i)=Rcum;
    Vout=[Vout;Vclass(2:end,:)];
    Snout=[Snout;Snclass(2:end,:)];

%         if hospInc==0
%             Rt(i)=heComputeEigs(pr,beta,D,NNfeed,ntot,Sclass(end,:)');
%             %Rt(i)=heComputeEigs(pr,beta,D,NNfeed,ntot,y0(1:lx+4));
%         end
    %Rt(i)=heComputeEigs(pr,beta,D,NNfeed,ntot,Sclass(end,:)');

    t0=tend;
    if i<lt-1

        Xh2w=NNvec(1:lx,i+1)-NNvec(1:lx,i);%Addition to each wp next intervention step

        Xw2h=-Xh2w; 
        Xw2h(Xw2h<0)=0;
        Xw2h=Xw2h./NNvec(1:lx,i);
        Xw2h(NNvec(1:lx,i)==0)=0;

        if NNvec(lx+adInd,i)>0%when would this not be the case?
            Xh2w(Xh2w<0)=0;
            Xh2w=Xh2w/NNvec(lx+adInd,i);
        else
            Xh2w=0;
        end

        %Move all infection statuses:

        y0=reshape(y0,[ntot,nc]);%IC
        y0w2h=y0(1:lx,:).*repmat(Xw2h,1,nc);%IC%number of people to be put at home (+)
        y0w2h=[-y0w2h;sum(y0w2h,1)];

        y0h2w=y0(lx+adInd,:);
        y0h2w=kron(y0h2w,Xh2w);
        y0h2w=[y0h2w;-sum(y0h2w,1)];
        y0([1:lx,lx+adInd],:)=y0([1:lx,lx+adInd],:)+y0w2h+y0h2w;

        y0=reshape(y0,ntot*nc,1);

    end

    if p2.sw==1
        Rt(end)=[];
        break;
    end

end

%% SWITCHING:

if p2.sw==1

i=2;
j=2;

tvec=tvec(~isnan(tvec));
tend=tvec(end);

while toutAll(end)<tend 

    NNfeed=NNvec(:,i);
    NNfeed(NNfeed==0)=1;
    D=Dvec(:,:,i);

    %Vaccination Rollout by Sector
    NNnext=NNvec(:,i);
    NNnext(lx+[1,2])=1;
    NNnext([1:lx,lx+3])=NNnext([1:lx,lx+3])/sumWorkingAge;
    NNnext(end)=1;
    p2.ratep1=NNnext.*[repmat(p2.aratep1(3),lx,1);p2.aratep1];    
    p2.ratep2=NNnext.*[repmat(p2.aratep2(3),lx,1);p2.aratep2];
    p2.ratep3=NNnext.*[repmat(p2.aratep3(3),lx,1);p2.aratep3];
    p2.ratep4=NNnext.*[repmat(p2.aratep4(3),lx,1);p2.aratep4];

    [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0,Sv1class,Vclass,Snclass]=integr8(data,NNfeed,D,i,t0,tend,dis,seedvec,y0,p2);

    toutAll=[toutAll;tout(2:end)];
    Sout=[Sout;Sclass(2:end,:)];
    Sv1out=[Sv1out;Sv1class(2:end,:)];
    Iout=[Iout;Itot(2:end)];
    Hout=[Hout;Hclass(2:end,:)];
    Dout=[Dout;Dclass(2:end,:)];
    DEout(:,j)=DEcum;
    %Rout(:,j)=Rcum;
    Vout=[Vout;Vclass(2:end,:)];
    Snout=[Snout;Snclass(2:end,:)];

%         if hospInc==0
%             Rt(j)=heComputeEigs(pr,beta,D,NNfeed,ntot,Sclass(end,:)');
%         end

    if toutAll(end)<tend

        tvec=[tvec(1:end-1),toutAll(end),tend];

        t0=toutAll(end);

        if i==2
        Xh2w=NNvec(1:lx,i+1)-NNvec(1:lx,i);%Addition to each wp next intervention step
        elseif i==3
        Xh2w=NNvec(1:lx,i-1)-NNvec(1:lx,i);%Addition to each wp next intervention step    
        end

        Xw2h=-Xh2w; 
        Xw2h(Xw2h<0)=0;
        Xw2h=Xw2h./NNvec(1:lx,i);
        Xw2h(NNvec(1:lx,i)==0)=0;

        if NNvec(lx+adInd,i)>0%when would this not be the case?
            Xh2w(Xh2w<0)=0;
            Xh2w=Xh2w/NNvec(lx+adInd,i);
        else
            Xh2w=0;
        end

        %Move all infection statuses:

        y0=reshape(y0,[ntot,nc]);%IC
        y0w2h=y0(1:lx,:).*repmat(Xw2h,1,nc);%IC%number of people to be put at home (+)
        y0w2h=[-y0w2h;sum(y0w2h,1)];

        y0h2w=y0(lx+adInd,:);
        y0h2w=kron(y0h2w,Xh2w);
        y0h2w=[y0h2w;-sum(y0h2w,1)];
        y0([1:lx,lx+adInd],:)=y0([1:lx,lx+adInd],:)+y0w2h+y0h2w;

        y0=reshape(y0,ntot*nc,1);

        if i==2
            i=i+1; 
        elseif i==3
            i=i-1; 
        end
        j=j+1;

    end   

end

end
    
%% OUTPUTS:  
    
f=[toutAll,...
   sum(Sout,2),...
   Iout,...
   sum(Hout,2),...
   sum(Dout,2),...
   sum(Vout(:,lx+1),2),...
   sum(Vout(:,lx+2),2),...
   sum(Vout(:,[1:lx,lx+3]),2),...
   sum(Vout(:,lx+4),2),...
   sum(Dout(:,lx+1),2),...
   sum(Dout(:,lx+2),2),...
   sum(Dout(:,[1:lx,lx+3]),2),...
   sum(Dout(:,lx+4),2),...
   sum(Snout,2)];%sum(DEout(end,:));
g=[max(sum(Hout(toutAll>=tvec(1),:),2)),...
   Rt(end)];

end

%%

function [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0new,Sv1class,Vclass,Snclass]=integr8(data,NN0,D,i,t0,tend,dis,seedvec,y0,p2)
%% CALL:

ntot=size(data.NNs,1);

fun=@(t,y)ODEs(data,NN0,D,t,dis,seedvec,y,p2);

if p2.sw==0||(p2.sw==1&&i==1)
    
    [tout,yout]=ode45(fun,[t0 tend],y0);
    %[tout,yout]=ode45(fun,(round(t0):1:tend),y0);
    
elseif p2.sw==1&&i==2
    
    %options=odeset('Events',@(t,y)changeinbehave(t,y,ntot,pr,p2));
    options=odeset('Events',@(t,y)lockdown(t,y,ntot,dis,p2));
    [tout,yout,te,ye,ie]=ode45(fun,[t0 tend],y0,options);
    %[tout,yout,te,ye,ie]=ode45(fun,(round(t0):1:tend),y0,options);
    
elseif p2.sw==1&&i==3
    
    %options=odeset('Events',@(t,y)changeinbehave(t,y,ntot,pr,p2));
    options=odeset('Events',@(t,y)reopen(t,y,ntot,dis,p2));
    [tout,yout,te,ye,ie]=ode45(fun,[t0 tend],y0,options);
    %[tout,yout,te,ye,ie]=ode45(fun,(round(t0):1:tend),y0,options);
    
end

%% EC:

Sclass=     yout(:,  0*ntot+1:1*ntot);
Sv1class=   yout(:,  8*ntot+1:10*ntot);
Itot=   sum(yout(:,  2*ntot+1:6*ntot),2)+sum(yout(:,  11*ntot+1:15*ntot),2);
Hclass=     yout(:,  6*ntot+1:7*ntot)       +yout(:,  15*ntot+1:16*ntot);
Rcum=       yout(end,7*ntot+1:8*ntot)       +yout(end,16*ntot+1:17*ntot);
Dclass=     yout(:,  17*ntot+1:18*ntot);
DEcum=      yout(end,17*ntot+1:18*ntot);
Vclass=     yout(:,18*ntot+1:19*ntot);
Snclass=    yout(:,  19*ntot+1:20*ntot);
y0new=      yout(end,:)';

end

%%

function [f,g]=ODEs(data,NN0,D,t,dis,seedvec,y,p2)

ntot=size(data.NNs,1);

%% IC:

S=      y(0*ntot+1:1*ntot);
E=      y(1*ntot+1:2*ntot);
Ina=    y(2*ntot+1:3*ntot);
Isa=    y(3*ntot+1:4*ntot);
Ins=    y(4*ntot+1:5*ntot);
Iss=    y(5*ntot+1:6*ntot);
H=      y(6*ntot+1:7*ntot);
R=      y(7*ntot+1:8*ntot);

Shv1=   y(8*ntot+1:9*ntot);
Sv1=    y(9*ntot+1:10*ntot);
Ev1=    y(10*ntot+1:11*ntot);
Inav1=  y(11*ntot+1:12*ntot);
Isav1=  y(12*ntot+1:13*ntot);
Insv1=  y(13*ntot+1:14*ntot);
Issv1=  y(14*ntot+1:15*ntot);
Hv1=    y(15*ntot+1:16*ntot);
Rv1=    y(16*ntot+1:17*ntot);

DE=     y(17*ntot+1:18*ntot);
V=      y(18*ntot+1:19*ntot);
Sn=     y(19*ntot+1:20*ntot);

%% TIME-DEPENDENT DISEASE PARAMETERS:

%Amplitude
amp = (Sn+(1-dis.heff).*(S-Sn))./S;

%Probabilities
ph = amp.*dis.ph;

%Calculations
Ts = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);

sig1 = dis.sig1;
sig2 = dis.sig2;
g1   = dis.g1;
g2   = (1-ph)./Ts;
g3   = dis.g3;
h    = ph./Ts;
mu   = dis.mu;
nu   = dis.nu;

%Transmission
red  = dis.red;
beta = dis.beta;

%Vaccination
hrv1  = dis.hrv1;    
scv1  = dis.scv1;  
g2_v1 = dis.g2_v1;
g3_v1 = dis.g3_v1;
h_v1  = dis.h_v1;
mu_v1 = dis.mu_v1;
trv1  = dis.trv1;
nuv1  = dis.nuv1;

%Preparedness
Tm      = p2.Tm;
p3      = p2.p3;
p4      = p2.p4;
%betamod = p2.betamod;
Hmax    = p2.Hmax;
g3_oc   = p2.g3_oc;
g3_ocv1 = p2.g3_ocv1;
mu_oc   = p2.mu_oc;
mu_ocv1 = p2.mu_ocv1;

startp1 = p2.startp1;
startp2 = p2.startp2;
startp3 = p2.startp3;
startp4 = p2.startp4;
pend    = p2.end;

ratep1 = p2.ratep1;
ratep2 = p2.ratep2;
ratep3 = p2.ratep3;
ratep4 = p2.ratep4;

%% FOI:

phi=1+data.amp*cos((t-32-data.phi)/(365/2*pi));

I=(red*Ina+Ins)  +(1-trv1)*(red*Inav1+Insv1);%Only non-self-isolating compartments

foi=phi*beta*(D*(I./NN0));

seed=phi*(seedvec./NN0);

%% HOSPITAL OCCUPANCY:

occ=max(1,sum(H+Hv1));

%% SELF-ISOLATION:

if t<Tm
    
    p3=0;
    p4=0;
    %p5=0;
    
else
    
    p3=p3;
    p4=p4;
    %p5=p5;
    
end

%% VACCINATION:

% %Uptake
% uptake=uptake-V./NNage;
% uptake(uptake<0)=0;
% uptake(uptake>0)=1;
% uptake=[repmat(uptake(3),numSectors,1);uptake];

nonVax=NN0-V;
%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised individuals
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%

if t>=pend
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
% elseif t>=p2.startp5
%     v1rates=p2.ratep5.*S./nonVax;
%     v1rater=p2.ratep5.*R./nonVax;
%     Vdot=   p2.ratep5;
      
elseif t>=startp4
    v1rates=ratep4.*S./nonVax;
    v1rater=ratep4.*R./nonVax;
    Vdot=   ratep4;
    
elseif t>=startp3
    v1rates=ratep3.*S./nonVax;
    v1rater=ratep3.*R./nonVax;
    Vdot=   ratep3;
    
elseif t>=startp2
    v1rates=ratep2.*S./nonVax;
    v1rater=ratep2.*R./nonVax;
    Vdot=   ratep2;
    
elseif t>=startp1
    v1rates=ratep1.*S./nonVax;
    v1rater=ratep1.*R./nonVax;
    Vdot=   ratep1;
    
else
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
end

%v1rate=0.2*p2.ratep3.*S./nonVax;
%v2rate=0.8*p2.ratep3.*S./nonVax;

Sndot=      -Sn.*(foi+seed)          -v1rates.*Sn./S;

%% EQUATIONS:

Sdot=       -S.*(foi+seed)  +nu.*R   -v1rates                +nuv1.*Sv1;
Shv1dot=    v1rates     -hrv1*Shv1   -Shv1.*foi;
Sv1dot=                  hrv1*Shv1   -Sv1.*(1-scv1).*foi  -nuv1.*Sv1;  %+nu.*Rv1; 

Edot=        S.*(foi+seed)  +Shv1.*foi  -(sig1+sig2).*E;
Ev1dot=                                  Sv1.*(1-scv1).*foi  -(sig1+sig2).*Ev1;

Inadot=     (1-p3)                   *sig1   .*E         -g1.*Ina;
Isadot=     p3                       *sig1   .*E         -g1.*Isa;
Insdot=     (1-p4)                   *sig2   .*E         -(g2+h)               .*Ins;
Issdot=     p4                       *sig2   .*E         -(g2+h)         .*Iss;
Inav1dot=   (1-p3)        *sig1   .*Ev1       -g1          .*Inav1;
Isav1dot=   p3            *sig1   .*Ev1       -g1         .*Isav1;
Insv1dot=   (1-p4)      *sig2   .*Ev1       -(g2_v1+h_v1)         .*Insv1;
Issv1dot=   p4          *sig2   .*Ev1       -(g2_v1+h_v1)   .*Issv1;

Hdot=       h.*(Ins+Iss)    -(g3+mu).*(min(occ,Hmax)/occ).*H   -(g3_oc+mu_oc).*(max(0,occ-Hmax)/occ).*H;
Hv1dot=     h_v1.*(Insv1+Issv1)  -(g3_v1+mu_v1).*(min(occ,Hmax)/occ).*Hv1   -(g3_ocv1+mu_ocv1).*(max(0,occ-Hmax)/occ).*Hv1;

Rdot= g1.*(Ina+Isa)+g2.*(Ins+Iss) + g3.*(min(occ,Hmax)/occ).*H  + g3_oc.*(max(0,occ-Hmax)/occ).*H    -nu.*R   -v1rater;
Rv1dot= g1.*(Inav1+Isav1)+g2_v1.*(Insv1+Issv1)+g3_v1.*(min(occ,Hmax)/occ).*Hv1 +g3_ocv1.*(max(0,occ-Hmax)/occ).*Hv1+v1rater;%-nu.*Rv1;

DEdot= mu.*(min(occ,Hmax)/occ).*H+mu_oc .*(max(0,occ-Hmax)/occ).*H+mu_v1.*(min(occ,Hmax)/occ).*Hv1  +mu_ocv1.*(max(0,occ-Hmax)/occ).*Hv1;     

%% OUTPUT:

f= [Sdot;Edot;...
    Inadot;Isadot;Insdot;Issdot;...
    Hdot;Rdot;...
    Shv1dot;Sv1dot;Ev1dot;...
    Inav1dot;Isav1dot;Insv1dot;Issv1dot;...
    Hv1dot;Rv1dot;...
    DEdot;Vdot;Sndot];
% f= [Sdot;Edot;...
%     Iadot;Ipdot;...
%     Inmdot;Ismdot;...
%     Insdot;Issdot;...
%     Qmdot;Qsdot;...
%     Hdot;DEdot;Rdot];

g=h.*(Ins+Iss)+h_v1.*(Insv1+Issv1);%Hin

end

%%

function [value,isterminal,direction]=changeinbehave(t,y,ntot,pr,p2)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = min(occ-p2.Hmax,0)+min(t-p2.Tm-0.1,0);
    direction  = 0;
    isterminal = ones(size(value));
    
end

function [value,isterminal,direction]=lockdown(t,y,ntot,dis,p2)
    
    S=      y(0*ntot+1:1*ntot);
    Ins=    y(4*ntot+1:5*ntot);
    Iss=    y(5*ntot+1:6*ntot);
    H=      y(6*ntot+1:7*ntot);
    Insv1=  y(13*ntot+1:14*ntot);
    Issv1=  y(14*ntot+1:15*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    Sn=     y(19*ntot+1:20*ntot);
    occ=    max(1,sum(H+Hv1));
    
    g3   = dis.g3;
    amp  = (Sn+(1-dis.heff).*(S-Sn))./S;
    ph   = amp.*dis.ph;
    Ts   = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    h    = ph./Ts;
    mu   = dis.mu;
    
    g3_v1 = dis.g3_v1;
    h_v1  = dis.h_v1;
    mu_v1 = dis.mu_v1;
    
    Hmax    = p2.Hmax;
    g3_oc   = p2.g3_oc;
    g3_ocv1 = p2.g3_ocv1;
    mu_oc   = p2.mu_oc;
    mu_ocv1 = p2.mu_ocv1;
    
    Hdot=       h.*(Ins+Iss)    -(g3+mu).*(min(occ,Hmax)/occ).*H   -(g3_oc+mu_oc).*(max(0,occ-Hmax)/occ).*H;
    Hv1dot=     h_v1.*(Insv1+Issv1)  -(g3_v1+mu_v1).*(min(occ,Hmax)/occ).*Hv1   -(g3_ocv1+mu_ocv1).*(max(0,occ-Hmax)/occ).*Hv1;
    occdot= sum(Hdot+Hv1dot);
     
    r=occdot/occ;
    Tcap=t+log(p2.Hmax/occ)/r;
    Tcap=Tcap-3;

    value      = [min(t-Tcap,0)+min(t-p2.Tm-0.1,0),min(occ-0.95*p2.Hmax,0)+min(t-p2.Tm-0.1,0)];
    direction  = [1,1];
    if r>0.025
        isterminal(1) = 1;
    else
        isterminal(1) = 0;
    end
    isterminal(2) = 1;
    
end

function [value,isterminal,direction]=reopen(t,y,ntot,dis,p2)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = occ-p2.thl;
    direction  = -1;
    isterminal = ones(size(value));
    
end
function [f,g,tvec]=heSimCovid19vax(pr,vx,beta,tvec,Dvec,n,ntot,NNvec,phi1,phi2,seedvec,S0,tau,plotTau)
%% PARAMETERS:

adInd=3;
lx=ntot-4;
NNbar=NNvec(:,1);
sumWorkingAge=sum(NNbar([1:lx,lx+3]));

%{
%Feed in to function - from prep
NNages=NNvec(:,1);
NNages=reshape(NNages,n,na);
NNages=sum(NNages,1);
%}
%DEout,Rout
%sigma,omega,gamma,hvec,muvec,pvec,qvec,n,ntot,na,NN,NNbar,NNrep,minNind,maxNind,maxN,Kbar,K1,Cbar,betaS,betaI,betaD,beta3

nc=19;
solvetype=2;
hospInc=0;

%% DETERMINISTIC OR STOCHASTIC SOLVER:

if solvetype==2
    
    lt=length(tvec);
    zn=zeros(ntot,1);
    
    t0=tvec(1);
    y0=[S0;repmat(zn,6,1);NNbar-S0;repmat(zn,nc-8,1)];
    
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
        vx.ratep1=NNnext.*[repmat(vx.aratep1(3),lx,1);vx.aratep1];    
        vx.ratep2=NNnext.*[repmat(vx.aratep2(3),lx,1);vx.aratep2];
        vx.ratep3=NNnext.*[repmat(vx.aratep3(3),lx,1);vx.aratep3];
        vx.ratep4=NNnext.*[repmat(vx.aratep4(3),lx,1);vx.aratep4];
        vx.ratep5=NNnext.*[repmat(vx.aratep5(3),lx,1);vx.aratep5];
        
        %pr.z=0;
        %{
        if i>=2%3
            pr.q1=pr.qnew;%Won't turn off****
            pr.q2=pr.qnew;
            %pr.z=pr.znew;
            pr.odds=pr.oddsnew;
        end
        %}
        %topen=1;%tvec(3);
        
        if hospInc>=1
            Rt(i)=heComputeEigs(pr,beta,D,NNfeed,ntot,y0(1:lx+4));
        end
        
        [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0,Hnew,Sv1class,Vclass]=integr8(pr,vx,beta,ntot,NNfeed,D,phi1,phi2,seedvec,t0,tend,y0,i,hospInc);
        
        toutAll=[toutAll;tout(2:end)];
        Sout=[Sout;Sclass(2:end,:)];
        Sv1out=[Sv1out;Sv1class(2:end,:)];
        Iout=[Iout;Itot(2:end)];
        Hout=[Hout;Hclass(2:end,:)];
        HnewAll=[HnewAll;Hnew(2:end)];
        Dout=[Dout;Dclass(2:end,:)];
        DEout(:,i)=DEcum;
        %Rout(:,i)=Rcum;
        Vout=[Vout;Vclass(2:end,:)];
               
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
        
        if pr.sw==1
            Rt(end)=[];
            break;
        end
        
    end
    
    %% SWITCHING:
    
    if pr.sw==1
    
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
        vx.ratep1=NNnext.*[repmat(vx.aratep1(3),lx,1);vx.aratep1];    
        vx.ratep2=NNnext.*[repmat(vx.aratep2(3),lx,1);vx.aratep2];
        vx.ratep3=NNnext.*[repmat(vx.aratep3(3),lx,1);vx.aratep3];
        vx.ratep4=NNnext.*[repmat(vx.aratep4(3),lx,1);vx.aratep4];
        
        [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0,Hnew,Sv1class,Vclass]=integr8(pr,vx,beta,ntot,NNfeed,D,phi1,phi2,seedvec,t0,tend,y0,i,hospInc);
                
        toutAll=[toutAll;tout(2:end)];
        Sout=[Sout;Sclass(2:end,:)];
        Sv1out=[Sv1out;Sv1class(2:end,:)];
        Iout=[Iout;Itot(2:end)];
        Hout=[Hout;Hclass(2:end,:)];
        HnewAll=[HnewAll;Hnew];
        Dout=[Dout;Dclass(2:end,:)];
        DEout(:,j)=DEcum;
        %Rout(:,j)=Rcum;
        Vout=[Vout;Vclass(2:end,:)];
               
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
    
elseif solvetype==1
    
    error('Final size calculations not possible')
    
elseif solvetype==3
    
    error('Code not written yet')
    %f=stochSim(y0,beta,gamma,n,ntot,NN,NN0,D,seed,phi1,phi2,tau,alpha);
    
end

%% OUTPUTS:
 
if tau==plotTau

    dodiff=1;
    plotEpi(toutAll,Sout,Hout,ntot,dodiff,tvec);%-tvec(2)+38

end    

if hospInc==0
    
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
       sum(Dout(:,lx+4),2)];%sum(DEout(end,:));
    %f=[toutAll(toutAll>0),sum(Hout(toutAll>0,:),2)];%For epi fit
    %Main constraints
    g=[max(sum(Hout(toutAll>=tvec(1),:),2)),...
       Rt(end)];
   %Rt calculated at end of each period
   
elseif hospInc==1%****
    
    f=[toutAll(toutAll>0),HnewAll(toutAll>0,:)];
    g=sum(sum(HnewAll(toutAll>tvec(3),:)));
    %g=Rt(3);
    
elseif hospInc==2
    
    f=[toutAll(toutAll>0),sum(Sout(toutAll>0,:),2),sum(Hout(toutAll>0,:),2)];%Occupancy
    g=sum(sum(HnewAll(toutAll>tvec(3),:)));%Incidence
    
end

end

%%

function [tout,Sclass,Hclass,Dclass,DEcum,Rcum,Itot,y0new,Hnew,Sv1class,Vclass]=integr8(pr,vx,beta,ntot,NN0,D,phi1,phi2,seedvec,t0,tend,y0,i,hospInc)
%% CALL:

fun=@(t,y)integr8covid(t,y,pr,vx,beta,ntot,NN0,D,phi1,phi2,seedvec);

if pr.sw==0||(pr.sw==1&&i==1)
    
    [tout,yout]=ode45(fun,[t0 tend],y0);
    %[tout,yout]=ode45(fun,(round(t0):1:tend),y0);
    
elseif pr.sw==1&&i==2
    
    options=odeset('Events',@(t,y)changeinbehave(t,y,ntot,pr,vx));
    %options=odeset('Events',@(t,y)lockdown(t,y,ntot,pr,vx));
    [tout,yout,te,ye,ie]=ode45(fun,[t0 tend],y0,options);
    %[tout,yout,te,ye,ie]=ode45(fun,(round(t0):1:tend),y0,options);
    
elseif pr.sw==1&&i==3
    
    options=odeset('Events',@(t,y)changeinbehave(t,y,ntot,pr,vx));
    %options=odeset('Events',@(t,y)reopen(t,y,ntot,pr,vx));
    [tout,yout,te,ye,ie]=ode45(fun,[t0 tend],y0,options);
    %[tout,yout,te,ye,ie]=ode45(fun,(round(t0):1:tend),y0,options);
    
end

% %% FOR FIT:   
% 
% if hospInc==1
%     Hnew=zeros(length(tout),ntot);
%     for i=1:length(tout)
%         [~,Hnewi]=fun(tout(i),yout(i,:)');
%         Hnew(i,:)=Hnewi';
%     end
% else
     Hnew=0;
% end

%% EC:

Sclass=     yout(:,  0*ntot+1:1*ntot);
Sv1class=   yout(:,  8*ntot+1:10*ntot);
Itot=   sum(yout(:,  2*ntot+1:6*ntot),2)+sum(yout(:,  11*ntot+1:15*ntot),2);
Hclass=     yout(:,  6*ntot+1:7*ntot)       +yout(:,  15*ntot+1:16*ntot);
Rcum=       yout(end,7*ntot+1:8*ntot)       +yout(end,16*ntot+1:17*ntot);
Dclass=     yout(:,  17*ntot+1:18*ntot);
DEcum=      yout(end,17*ntot+1:18*ntot);
Vclass=     yout(:,18*ntot+1:19*ntot);
y0new=      yout(end,:)';

end

%%

function [f,g]=integr8covid(t,y,pr,vx,betaIn,ntot,NN0,D,phi1,phi2,seedvec)
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

% Inm=  y(4*ntot+1:5*ntot);
% Ism=  y(5*ntot+1:6*ntot);
% Qm=   y(7*ntot+1:8*ntot);
% Qs=   y(8*ntot+1:9*ntot);
% H2=   y(11*ntot+1:12*ntot);

%% FOI:

phi=phi1;
%phi=phi1-phi2*cos(pi*t/180);%Seasonality****

beta=betaIn;
%{
if t>topen && t<topen+30
    beta=betaIn*(.25+.75*(t-topen)/30);
else
    beta=betaIn;
end
%}

I=(pr.red*Ina+Ins)  +(1-vx.trv1)*(pr.red*Inav1+Insv1);%Only non-self-isolating compartments

foi=phi*beta*(D*(I./NN0));

seed=phi*(seedvec./NN0);

%% HOSPITAL OCCUPANCY:

occ=max(1,sum(H+Hv1));

%% SELF-ISOLATION:

if t<pr.Tm
    
    p3=0;
    p4=0;
    %p5=0;
    
else
    
    p3=pr.p3;
    p4=pr.p4;
    %p5=pr.p5;
    
end

%% VACCINATION:

% %Uptake
% uptake=vx.uptake-V./vx.NNage;
% uptake(uptake<0)=0;
% uptake(uptake>0)=1;
% uptake=[repmat(uptake(3),numSectors,1);uptake];

nonVax=NN0-V;
%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised individuals
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%

if t>=vx.end
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
elseif t>=vx.startp5
    v1rates=vx.ratep5.*S./nonVax;
    v1rater=vx.ratep5.*R./nonVax;
    Vdot=   vx.ratep5;
      
elseif t>=vx.startp4
    v1rates=vx.ratep4.*S./nonVax;
    v1rater=vx.ratep4.*R./nonVax;
    Vdot=   vx.ratep4;
    
elseif t>=vx.startp3
    v1rates=vx.ratep3.*S./nonVax;
    v1rater=vx.ratep3.*R./nonVax;
    Vdot=   vx.ratep3;
    
elseif t>=vx.startp2
    v1rates=vx.ratep2.*S./nonVax;
    v1rater=vx.ratep2.*R./nonVax;
    Vdot=   vx.ratep2;
    
elseif t>=vx.startp1
    v1rates=vx.ratep1.*S./nonVax;
    v1rater=vx.ratep1.*R./nonVax;
    Vdot=   vx.ratep1;
    
else
    v1rates=zeros(ntot,1);
    v1rater=zeros(ntot,1);
    Vdot=   zeros(ntot,1);
    
end

%v1rate=0.2*vx.ratep3.*S./nonVax;
%v2rate=0.8*vx.ratep3.*S./nonVax;

%% EQUATIONS:

Sdot=       -S.*(foi+seed)  +pr.nu.*R   -v1rates                +vx.nuv1.*Sv1;
Edot=        S.*(foi+seed)  +Shv1.*foi  -pr.sigma.*E;

Inadot=     (1-p3)      *(1-pr.p1)              *pr.sigma   .*E         -pr.g1*(1+pr.odds)          .*Ina;
Isadot=     p3          *(1-pr.p1)              *pr.sigma   .*E         -pr.g1*(1+pr.odds)          .*Isa;
Insdot=     (1-p4)      *pr.p1                  *pr.sigma   .*E         -(pr.g2+pr.h)               .*Ins;
Issdot=     p4          *pr.p1                  *pr.sigma   .*E         -(pr.g2+pr.h+pr.q2)         .*Iss;

Hdot=       pr.h.*(Ins+Iss)                                             -(pr.g3+pr.mu).*(min(occ,pr.Hmax)/occ).*H           -(pr.g3_oc+pr.mu_oc).*(max(0,occ-pr.Hmax)/occ).*H;
Rdot=       pr.g1.*(1+pr.odds).*(Ina+Isa)+pr.g2.*(Ins+Iss)              + pr.g3.*(min(occ,pr.Hmax)/occ).*H                  + pr.g3_oc.*(max(0,occ-pr.Hmax)/occ).*H         -pr.nu.*R   -v1rater;

Shv1dot=    v1rates     -vx.hrv1*Shv1   -Shv1.*foi;
Sv1dot=                  vx.hrv1*Shv1   -Sv1.*(1-vx.scv1).*foi  -vx.nuv1.*Sv1;  %+pr.nu.*Rv1;   
Ev1dot=                                  Sv1.*(1-vx.scv1).*foi  -pr.sigma.*Ev1;

Inav1dot=   (1-p3)      *(1-pr.p1*(1-vx.p1v1))  *pr.sigma   .*Ev1       -pr.g1*(1+pr.odds)          .*Inav1;
Isav1dot=   p3          *(1-pr.p1*(1-vx.p1v1))  *pr.sigma   .*Ev1       -pr.g1*(1+pr.odds)          .*Isav1;
Insv1dot=   (1-p4)      *(1-vx.p1v1)*pr.p1      *pr.sigma   .*Ev1       -(vx.g2_v1+vx.h_v1)         .*Insv1;
Issv1dot=   p4          *(1-vx.p1v1)*pr.p1      *pr.sigma   .*Ev1       -(vx.g2_v1+vx.h_v1+pr.q2)   .*Issv1;

Hv1dot=     vx.h_v1.*(Insv1+Issv1)                                      -(vx.g3_v1+vx.mu_v1).*(min(occ,pr.Hmax)/occ).*Hv1   -(vx.g3_ocv1+vx.mu_ocv1).*(max(0,occ-pr.Hmax)/occ).*Hv1;
Rv1dot=     pr.g1.*(1+pr.odds).*(Inav1+Isav1)   +vx.g2_v1.*(Insv1+Issv1)+ vx.g3_v1.*(min(occ,pr.Hmax)/occ).*Hv1             + vx.g3_ocv1.*(max(0,occ-pr.Hmax)/occ).*Hv1     +v1rater;%-pr.nu.*Rv1;

DEdot=      pr.mu.*(min(occ,pr.Hmax)/occ).*H    +pr.mu_oc .*(max(0,occ-pr.Hmax)/occ).*H     +vx.mu_v1.*(min(occ,pr.Hmax)/occ).*Hv1  +vx.mu_ocv1.*(max(0,occ-pr.Hmax)/occ).*Hv1;     

% Inmdot= (1-p5)      *pr.p1      *(1-pr.p2)  *pr.sigma   .*E           -pr.gX                  .*Inm;
% Ismdot= p5          *pr.p1      *(1-pr.p2)  *pr.sigma   .*E           -(pr.gX+pr.q1)          .*Ism;
% Ipdot=pr.p1*pr.sigma*E-(1+pr.odds)*pr.omega*Ip;%2**
% Qmdot=  pr.g1*pr.odds*Ia+pr.q1*Ism-pr.g4.*Qm;
% Qsdot=  pr.q2*Iss-pr.h.*Qs;%no recovery rate here???

%% OUTPUT:

f= [Sdot;Edot;...
    Inadot;Isadot;Insdot;Issdot;...
    Hdot;Rdot;...
    Shv1dot;Sv1dot;Ev1dot;...
    Inav1dot;Isav1dot;Insv1dot;Issv1dot;...
    Hv1dot;Rv1dot;...
    DEdot;Vdot];
% f= [Sdot;Edot;...
%     Iadot;Ipdot;...
%     Inmdot;Ismdot;...
%     Insdot;Issdot;...
%     Qmdot;Qsdot;...
%     Hdot;DEdot;Rdot];

g=pr.h.*(Ins+Iss)+vx.h_v1.*(Insv1+Issv1);%Hin

end

%{
Ip=y(3*ntot+1:4*ntot);
Inm=y(4*ntot+1:5*ntot);
Ism=y(5*ntot+1:6*ntot);
Ins=y(6*ntot+1:7*ntot);
Iss=y(7*ntot+1:8*ntot);
Qm=y(8*ntot+1:9*ntot);
Qs=y(9*ntot+1:10*ntot);
H=y(10*ntot+1:11*ntot);
%H2=y(11*ntot+1:12*ntot);
I=2/3*Ia+2/3*Ip+Inm+Ism+Ins+Iss;%All infectious
%DE=y(10*ntot+1:11*ntot);
%R=y(11*ntot+1:end);
%}
%{
%HE model:
Iadot=(1-pr.p1)*pr.sigma*E-(1+pr.odds)*pr.g1*Ia;%2**
Ipdot=pr.p1*pr.sigma*E-(1+pr.odds)*pr.omega*Ip;%2**
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-pr.g2*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-pr.h*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.g1*pr.odds*Ia+(1-pr.p2)*pr.omega*pr.odds.*Ip+pr.q1*Ism-pr.g4*Qm;
Qsdot=pr.p2*pr.omega*pr.odds.*Ip+pr.q2*Iss-pr.h*Qs;
%}
%{
Iadot=(1-pr.p1)*pr.sigma*E-pr.g1*Ia;
Ipdot=pr.p1*pr.sigma*E-pr.omega*Ip;
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-pr.g2*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-pr.h*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.q1*Ism-pr.g4*Qm;
Qsdot=pr.q2*Iss-pr.h*Qs;
%}
%{
%With quarantine of all infs:
Iadot=(1-pr.p1)*pr.sigma*E-(pr.g1+pr.q1*pr.z)*Ia;
Ipdot=pr.p1*pr.sigma*E-(pr.omega+pr.q1*pr.z)*Ip;%q1=q2****
Inmdot=(1-pr.p2)*(1-pr.p3)*pr.omega.*Ip-(pr.g2+pr.q1*pr.z)*Inm;
Ismdot=(1-pr.p2)*pr.p3*pr.omega.*Ip-(pr.g2+pr.q1)*Ism;
Insdot=pr.p2*(1-pr.p4)*pr.omega.*Ip-(pr.h+pr.q2*pr.z)*Ins;%*(1-h).*Ins;
Issdot=pr.p2*pr.p4*pr.omega.*Ip-(pr.h+pr.q2)*Iss;%(1-h).*Iss;
Qmdot=pr.q1*(Ia*pr.z+(1-pr.p2)*pr.z.*Ip+Inm*pr.z+Ism)-pr.g4*Qm;
Qsdot=pr.q2*(pr.p2.*Ip*pr.z+Ins*pr.z+Iss)-pr.h*Qs;
%}

%%

%Stochastic variant - needs update to C19 flowchart
function f=stochSim(y,beta,gamma,n,ntot,NN,N0,D,seed,phi1,phi2,tau,alpha)
%Still flu-like/SIR****
%Feed in mu if required
factor=6;
tend=360*factor; beta=beta/factor; gamma=gamma/factor;
Vec=zeros(ntot,tend);
%
S=y(1:ntot);
I=y(ntot+1:2*ntot);
R=y(2*ntot+1:end);
i=1;
threshold=30;%Number of time-steps for necessary simulation/seed
while i<tend && (i<threshold || sum(I)>0)%At least 30 time-steps
phi=1;%phi1-phi2*cos(pi*i*f1/180);%Just seed here
Sout=1-exp(-phi*(beta*(D*(I./N0).^alpha)+seed*heaviside(threshold-i)));%+mu*R;.^alpha
Sout(Sout>1)=1;
Sout=binornd(S,Sout); Sout(S==0)=0;
S=S-Sout; I=I+Sout;
%
Iout=1-exp(-gamma);
Iout=binornd(I,Iout);
I=I-Iout; R=R+Iout;
Vec(:,i)=I;
i=i+1;
end
f=R;
if sum(isnan(R))>0
    fuck=1;
    print('Somenting is NaN')
end
end

%%

function [value,isterminal,direction]=changeinbehave(t,y,ntot,pr,vx)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = min(occ-pr.Hmax,0)+min(t-pr.Tm-0.1,0);
    direction  = 0;
    isterminal = ones(size(value));
    
end

function [value,isterminal,direction]=lockdown(t,y,ntot,pr,vx)

    Ins=    y(4*ntot+1:5*ntot);
    Iss=    y(5*ntot+1:6*ntot);
    H=      y(6*ntot+1:7*ntot);
    Insv1=  y(13*ntot+1:14*ntot);
    Issv1=  y(14*ntot+1:15*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1)); 
    Hdot=   pr.h.*(Ins+Iss)         -(pr.g3+pr.mu).*(min(occ,pr.Hmax)/occ).*H           -(pr.g3_oc+pr.mu_oc).*(max(0,occ-pr.Hmax)/occ).*H;
    Hv1dot= vx.h_v1.*(Insv1+Issv1)  -(vx.g3_v1+vx.mu_v1).*(min(occ,pr.Hmax)/occ).*Hv1   -(vx.g3_ocv1+vx.mu_ocv1).*(max(0,occ-pr.Hmax)/occ).*Hv1;
    occdot= sum(Hdot+Hv1dot);
     
    r=occdot/occ;
    Tcap=t+log(pr.Hmax/occ)/r;
    Tcap=Tcap-7;

    value      = [min(t-Tcap,0)+min(t-pr.Tm-0.1,0),min(occ-0.95*pr.Hmax,0)+min(t-pr.Tm-0.1,0)];
    direction  = [1,1];
    if r>0.025
        isterminal(1) = 1;
    else
        isterminal(1) = 0;
    end
    isterminal(2) = 1;
    
end

function [value,isterminal,direction]=reopen(t,y,ntot,pr,vx)

    H=      y(6*ntot+1:7*ntot);
    Hv1=    y(15*ntot+1:16*ntot);
    occ=    max(1,sum(H+Hv1));

    value      = occ-pr.thl;
    direction  = -1;
    isterminal = ones(size(value));
    
end

%%

function f=plotEpi(tout,Y,H,n,dodiff,tvec)
yvar='Inc./hosp.';%'Susceptibles'; 'Hospitalisations';
solvetype=2;
tend=tvec(end);%720;%For plot only
na=size(Y,2)/n;
cmap=lines(7);
if solvetype==2
    if dodiff==1
        Y=-diff(Y,1);
        tdiff=diff(tout,1);
        Y=Y./repmat(tdiff,1,n*na);%repmat - older version?
        tout(1)=[];
        H(1,:)=[];
    end
    figure
    fs=10; lw=2;
    if na==4
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:end),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2),sum(H(:,3*n+1:end),2)];
        %Yall=Y(:,1:n)+Y(:,n+1:2*n)+Y(:,2*n+1:3*n)+Y(:,3*n+1:end);
    elseif na==5
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2),sum(Y(:,3*n+1:4*n),2),sum(Y(:,4*n+1:end),2)];
    elseif na==3
        Yall=[sum(Y(:,1:n),2),sum(Y(:,n+1:2*n),2),sum(Y(:,2*n+1:3*n),2)];
        Hall=[sum(H(:,1:n),2),sum(H(:,n+1:2*n),2),sum(H(:,2*n+1:3*n),2)];
    elseif na==1
        Yall=sum(Y,2);
        Hall=sum(H,2);
    else
        error('Number of age groups not recognised for plotting')
    end
    %maxY=max(max(Yall));
    maxY=max(max(max(Yall)),max(max(Hall)));
    %
    %Unlogged plots:
    h=zeros(1,na);
    hold on
    for i=2:length(tvec)-1
        plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
    end
    for i=1:na
        %h(i)=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(i,:));
        h1=plot(tout,Yall(:,i),'linewidth',lw,'color',cmap(1,:));%(i,:)
        h2=plot(tout,Hall(:,i),'-','linewidth',lw,'color',cmap(2,:));%'--', (i,:)
    end
    %plot(tout,Yall,'linewidth',.5,'color',cmap(2,:));
    %plot(tout,sum(Yall,2),'linewidth',lw,'color','k');
    %}
    %{
    %Logged plots:
    hold on
    semilogy(tout,Yall);
    %}
    lt=length(tvec);
    points=[0,tvec(2:end)]+10;
    pointsy=.93*maxY;
    txt={'PRE','LD','1','2','3','4','5','6'};
    for i=1:lt-1
        text(points(i),pointsy,txt{i},'fontsize',20)
    end
    
    xlabel('Time','FontSize',fs);
    ylabel('Population','FontSize',fs);%yvar
    set(gca,'FontSize',fs);
    axis([0,tend,0,maxY])
    
    xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
    xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
    
    legend([h1,h2],'Inc.','HO','location','W')
    %
    if na==4
        legend(h,{'0-4','5-19','20-64','65+'},'location','NE')
    elseif na==5
        legend(h,{'0-4','5-17','18-49','50-64','65+'},'location','NE')
    elseif na==3
        legend(h,{'0-15','16-64','65+'},'location','NE')
    end
    %}
    grid on
    grid minor
    box on
    hold off
%}
end
end
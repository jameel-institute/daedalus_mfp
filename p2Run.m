function [data,f,g]=p2Run(data,dis,p2);

ln    = length(data.NNs);
lx    = length(data.obj);
int   = length(data.xoptim)/lx;

XitMat = reshape(data.xoptim,lx,int);
Dvec   = zeros(ln,ln,int);
for i = 1:int;
    Dvec(:,:,i) = p2MakeDs(data,XitMat(:,i),data.hw(i,:));
end
data.Dvec = Dvec;

data.Ev = dis.Ev.*XitMat(data.IntlInd,:);

[data,f,g] = p2SimVax(data,dis,data.NNs,XitMat,p2);%S0 = data.NNs, i.e. entirely susceptible population

end

%%

function [data,f,g]=p2SimVax(data,dis,S0,XitMat,p2)               
%% IC:

t0    = data.tvec(1);
ln    = length(data.NNs);
lx    = length(data.obj);
adInd = 3;
nc    = 20;
zn    = zeros(ln,1);
NN    = data.NNs;
y0    = [S0;repmat(zn,6,1);NN-S0;repmat(zn,nc-9,1);S0];

%initialise outputs
Tout       = t0;
Tsout      = dis.Ts';
Thout      = dis.Th';
phout      = dis.ph';   
Iout       = zn';
Isaout     = zn';
Isavout    = zn';
Issout     = zn';
Issvout    = zn';
Insout     = zn';
Hout       = zn';
Dout       = zn';
dDout      = zn';
Xout       = [];
hwout      = [];
aaout      = 0;
asout      = 0;
betamodout = 1;
Vout       = zn';

%% LOOP:

i = 1;

tend = data.tvec(end);

while Tout(end)<tend; 
    
    Xit = XitMat(:,i);    
    D   = data.Dvec(:,:,i);
    
    [tout,Ts,Th,ph,Iclass,Isaclass,Isavlass,Issclass,Issvlass,Insclass,Hclass,Dclass,dDclass,asc_a,asc_s,betamod,Vclass,y0,inext] = integr8(data,D,i,t0,tend,dis,y0,p2);
    
    if tout(end)<tend;
        data.tvec = [data.tvec(1:end-1),tout(end),tend];
        p2.Tres   = min(p2.Tres,tout(end));
        t0        = tout(end);
        i         = inext;
    end   

    Tout       = [Tout;tout(2:end)];  
    Tsout      = [Tsout;Ts(2:end,:)];  
    Thout      = [Thout;Th(2:end,:)];
    phout      = [phout;ph(2:end,:)];
    Iout       = [Iout;Iclass(2:end,:)];
    Isaout     = [Isaout;Isaclass(2:end,:)];
    Isavout    = [Isavout;Isavlass(2:end,:)];
    Issout     = [Issout;Issclass(2:end,:)];
    Issvout    = [Issvout;Issvlass(2:end,:)];
    Insout     = [Insout;Insclass(2:end,:)];
    Hout       = [Hout;Hclass(2:end,:)];
    Dout       = [Dout;Dclass(2:end,:)];
    dDout      = [dDout;dDclass(2:end,:)]; 
    X          = Xit'.*ones(length(tout),lx);
    Xout       = [Xout;X(1:end-1,:)];      
    hw         = data.hw(i,:).*ones(length(tout),lx);
    hwout      = [hwout;hw(1:end-1,:)];
    aaout      = [aaout;asc_a(2:end)];
    asout      = [asout;asc_s(2:end)];
    betamodout = [betamodout;betamod(2:end)];
    Vout       = [Vout;Vclass(2:end,:)];
    
end
    
%% OUTPUTS:  

Xout  = [Xout;Xout(end,:)];
hwout = [hwout;hwout(end,:)];
g     = [Tout,Xout,hwout,Isaout,Isavout,Issout,Issvout,Insout,Hout,Dout,Vout,betamodout,Tsout,Thout,phout];
f     = [Tout,...
         sum(Iout,2),...
         sum(Hout,2),...
         sum(Dout,2),...
         aaout,...
         asout,...
         betamodout,...  
         sum(Vout(:,lx+1),2),...
         sum(Vout(:,lx+2),2),...
         sum(Vout(:,[1:lx,lx+adInd]),2),...
         sum(Vout(:,lx+4),2),...
         sum(Dout(:,lx+1),2),...
         sum(Dout(:,lx+2),2),...
         sum(Dout(:,[1:lx,lx+adInd]),2),...
         sum(Dout(:,lx+4),2),...
         sum(dDout(:,lx+1),2),...
         sum(dDout(:,lx+2),2),...
         sum(dDout(:,[1:lx,lx+adInd]),2),...
         sum(dDout(:,lx+4),2)];
  
end

%%

function [tout,Ts,Th,ph,Iclass,Isaclass,Isavlass,Issclass,Issvlass,Insclass,Hclass,Dclass,dDclass,asc_a,asc_s,betamod,Vclass,y0new,inext] = integr8(data,D,i,t0,tend,dis,y0,p2)
%% CALL:

ln  = length(data.NNs);
fun = @(t,y)ODEs(data,D,i,t,dis,y,p2);

if strcmp(data.inp3,'No Closures');
    options = odeset('Events',@(t,y)unmitigated(t,y,data,dis,i,p2));
elseif strcmp(data.inp3,'School Closures');
    options = odeset('Events',@(t,y)reactive_closures(t,y,data,dis,i,p2));
elseif strcmp(data.inp3,'Economic Closures');
    options = odeset('Events',@(t,y)reactive_closures(t,y,data,dis,i,p2));
elseif strcmp(data.inp3,'Elimination');
	options = odeset('Events',@(t,y)elimination(t,y,data,dis,i,p2));
else
    error('Unknown Mitigation Strategy!');
end    

[tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);

y0new     = yout(end,:)'; 
if tout(end)<tend;
    inext = data.inext(ie(end));
else
    inext = NaN;
end

%% OUTPUT VARIABLES:

Iclass   = yout(:, 2*ln+1: 3*ln) + yout(:, 3*ln+1: 4*ln) + ...
           yout(:, 4*ln+1: 5*ln) + yout(:, 5*ln+1: 6*ln) + ...
           yout(:,11*ln+1:12*ln) + yout(:,12*ln+1:13*ln) + ...
           yout(:,13*ln+1:14*ln) + yout(:,14*ln+1:15*ln);
Isaclass = yout(:, 3*ln+1: 4*ln); 
Isavlass = yout(:,12*ln+1:13*ln);
Issclass = yout(:, 5*ln+1: 6*ln);
Issvlass = yout(:,14*ln+1:15*ln);
Insclass = yout(:, 4*ln+1: 5*ln) + yout(:,13*ln+1:14*ln);
Hclass   = yout(:, 6*ln+1: 7*ln) + yout(:,15*ln+1:16*ln);
Dclass   = yout(:,17*ln+1:18*ln);
Vclass   = yout(:,18*ln+1:19*ln);

%% TIME-DEPENDENT PARAMETERS:

occ   = max(1,sum(Hclass,2));
Hmax  = p2.Hmax;
th0   = max(1,1+1.87*((occ-Hmax)/(2*Hmax-Hmax)));

pd  = min(th0.*dis.pd',1);
Th  = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
mu  = pd./Th;
ddk = 10^5*sum(mu.*Hclass,2)/sum(data.Npop);

sd_fun = @(l,b,c,t,d) l + (1-l)*exp(-b*exp(-c.*t).*d);%here, t is time since response

if strcmp(data.inp3,'No Closures')||i==1;
    betamod = ones(size(occ));
elseif any(i==data.imand);
    betamod = min(sd_fun(p2.sdl,p2.sdb,p2.sdc,tout-p2.Tres,ddk), sd_fun(p2.sdl,p2.sdb,p2.sdc,14,2));
else
    betamod = sd_fun(p2.sdl,p2.sdb,p2.sdc,tout-p2.Tres,ddk);
end

S     = yout(:,0*ln+1:1*ln);
Sn    = yout(:,19*ln+1:20*ln);
Ins   = yout(:,4*ln+1:5*ln);
Iss   = yout(:,5*ln+1:6*ln);
Insv1 = yout(:,13*ln+1:14*ln);
Issv1 = yout(:,14*ln+1:15*ln);
H     = yout(:,6*ln+1:7*ln);
Hv1   = yout(:,15*ln+1:16*ln);

amp   = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
ph    = amp.*dis.ph';
Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g3    = (1-pd)./Th;
h     = ph./Ts;
h_v1  = dis.h_v1;

Hdot   = h.*Ins      +h.*Iss        -(g3+mu).*H;
Hv1dot = h_v1'.*Insv1 +h_v1'.*Issv1   -(g3+mu).*Hv1;
occdot = sum(Hdot+Hv1dot,2);
r      = occdot./occ;

dDclass = mu.*H + mu.*Hv1;     

trate  = p2.trate;
ps     = dis.ps;
siga   = dis.siga;
sigs   = dis.sigs;
E      = yout(:,1*ln+1:2*ln);
Ev1    = yout(:,10*ln+1:11*ln);
incid  = max(0,10^5*((siga+sigs)*sum(E+Ev1,2))/sum(data.Npop));

if i~=5;
    asc_s  = 1./(1+exp(0.7623+1.605*log10(incid)-1.416*log10(trate)));
    propCT = 1./(1+exp(2.159+1.697*log10(incid)));
    asc_a  = propCT.*asc_s + (1-propCT)*0;
    
    asc_a  = min(asc_a,(trate./incid).*(0*(1-propCT) + (1-ps)*propCT));
    asc_s  = min(asc_s,(trate./incid).*(1*(1-propCT) + ps*propCT));
    asc_a  = max(trate/10^5*(0*(1-propCT) + (1-ps)*propCT),asc_a);
    asc_s  = max(trate/10^5*(1*(1-propCT) + ps*propCT),asc_s);    

    asc_a  = asc_a.*(tout>p2.t_tit).*(tout<p2.end);    
    asc_s  = asc_s.*(tout>p2.t_tit).*(tout<p2.end);    

else
    asc_a  = zeros(size(tout));
    asc_s  = zeros(size(tout));

end

end

%%

function [f,g]=ODEs(data,D,i,t,dis,y,p2)

NN        = data.NNs;
ln        = length(NN);
NN(NN==0) = 1;

%% IC:

S     = y(0*ln+1:1*ln);
E     = y(1*ln+1:2*ln);
Ina   = y(2*ln+1:3*ln);
Isa   = y(3*ln+1:4*ln);
Ins   = y(4*ln+1:5*ln);
Iss   = y(5*ln+1:6*ln);
H     = y(6*ln+1:7*ln);
R     = y(7*ln+1:8*ln);

Shv1  = y(8*ln+1:9*ln);
Sv1   = y(9*ln+1:10*ln);
Ev1   = y(10*ln+1:11*ln);
Inav1 = y(11*ln+1:12*ln);
Isav1 = y(12*ln+1:13*ln);
Insv1 = y(13*ln+1:14*ln);
Issv1 = y(14*ln+1:15*ln);
Hv1   = y(15*ln+1:16*ln);
Rv1   = y(16*ln+1:17*ln);

DE    = y(17*ln+1:18*ln);
V     = y(18*ln+1:19*ln);
Sn    = y(19*ln+1:20*ln);

%% HOSPITAL OCCUPANCY:

occ   = max(1,sum(H+Hv1));
Hmax  = p2.Hmax;

%% TIME-DEPENDENT DISEASE PARAMETERS:

%Amplitudes
amp = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
th0 = max(1,1+1.87*((occ-Hmax)/(2*Hmax-Hmax)));

%Probabilities
ps = dis.ps;
ph = amp.*dis.ph;
pd = min(th0*dis.pd,1);

%Durations
Tlat = dis.Tlat;
Tinc = dis.Tinc;
Tay  = dis.Tay;
Ts   = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
Th   = ((1-pd).*dis.Threc)+(pd.*dis.Thd);

%Rates
siga = dis.siga;
sigs = dis.sigs;
g1   = dis.g1;
g2   = (1-ph)./Ts;
g3   = (1-pd)./Th;
h    = ph./Ts;
mu   = pd./Th;
nu   = dis.nu;

%Transmission
red  = dis.red;
beta = dis.beta;

%Vaccination
hrv1  = dis.hrv1;    
scv1  = dis.scv1;  
g2_v1 = dis.g2_v1;
h_v1  = dis.h_v1;
trv1  = dis.trv1;
nuv1  = dis.nuv1;

%Preparedness
trate   = p2.trate;
startp1 = p2.startp1;
startp2 = p2.startp2;
startp3 = p2.startp3;
startp4 = p2.startp4;
pend    = p2.end;
aratep1 = p2.aratep1;
aratep2 = p2.aratep2;
aratep3 = p2.aratep3;
aratep4 = p2.aratep4;

%% SELF-ISOLATION:

if t<p2.t_tit;   
    asc_a = 0;
    asc_s = 0;
    tm_a  = 1;
    tm_s  = 1;
    
elseif t<p2.end && i~=5;
    incid  = max(0,10^5*((siga+sigs)*sum(E+Ev1))/sum(data.Npop));
    asc_s  = 1/(1+exp(0.7623+1.605*log10(incid)-1.416*log10(trate)));
    propCT = 1/(1+exp(2.159+1.697*log10(incid)));
    asc_a  = propCT*asc_s + (1-propCT)*0;
    
    asc_a = min(asc_a,(trate/incid)*(0*(1-propCT) + (1-ps)*propCT));
    asc_s = min(asc_s,(trate/incid)*(1*(1-propCT) + ps*propCT));
    asc_a = max(trate/10^5*(0*(1-propCT) + (1-ps)*propCT),asc_a);
    asc_s = max(trate/10^5*(1*(1-propCT) + ps*propCT),asc_s);
    
    onsPCR_s = 11.3224-2.6260*log10(trate);
    onsPCR_c = onsPCR_s - 5.6304;
    Teff_c   = max(0,Tinc+onsPCR_c-Tlat);
    Teff_s   = max(0,Tinc+onsPCR_s-Tlat);
    mult_ac  = min(Teff_c,Tay)./Tay;
    mult_sc  = min(Teff_c,Ts)./Ts;
    mult_ss  = min(Teff_s,Ts)./Ts;
    
    tm_a = mult_ac;
    tm_s = mult_sc*propCT + mult_ss*(1-propCT);
    
else
    asc_a = 0;
    asc_s = 0;
    tm_a  = 1;
    tm_s  = 1;
        
end

sig1 = siga*(1-asc_a);
sig2 = sigs*(1-asc_s);
sig3 = siga*asc_a;
sig4 = sigs*asc_s;

%% VACCINATION:

% %Uptake
% uptake=uptake-V./NNage;
% uptake(uptake<0)=0;
% uptake(uptake>0)=1;
% uptake=[repmat(uptake(3),numSectors,1);uptake];

%nonVax=NN0-V;
nonVax=S+E+Ina+R;%only living, asymptomatic, non-isolating, unvaccinated individuals are eligible for vaccination
%S (and R) ./nonVax accounts for the (inefficient) administration of vaccines to exposed, infectious and hospitalised people
%nonVax approximates S+E+I+H+R (unvaccinated) and D (partially)
%S (or R) ./nonVax is approximately 1 (or 0) when prevalence is low but is closer to 0 (or 1) when prevalence is high
%nonVax is non-zero as long as uptake is less than 100%

if t>=pend
    v1rate_sn = zeros(ln,1);
    v1rate_sw = zeros(ln,1);
    v1rate_r  = zeros(ln,1);
    Vdot      = zeros(ln,1);
      
elseif t>=startp4
    v1rate_sn = aratep4.*Sn./nonVax;
    v1rate_sw = aratep4.*(S-Sn)./nonVax;
    v1rate_r  = aratep4.*R./nonVax;
    Vdot      = aratep4;
    
elseif t>=startp3
    v1rate_sn = aratep3.*Sn./nonVax;
    v1rate_sw = aratep3.*(S-Sn)./nonVax;
    v1rate_r  = aratep3.*R./nonVax;
    Vdot      = aratep3;
    
elseif t>=startp2
    v1rate_sn = aratep2.*Sn./nonVax;
    v1rate_sw = aratep2.*(S-Sn)./nonVax;
    v1rate_r  = aratep2.*R./nonVax;
    Vdot      = aratep2;
    
elseif t>=startp1
    v1rate_sn = aratep1.*Sn./nonVax;
    v1rate_sw = aratep1.*(S-Sn)./nonVax;
    v1rate_r  = aratep1.*R./nonVax;
    Vdot      = aratep1;
    
else
    v1rate_sn = zeros(ln,1);
    v1rate_sw = zeros(ln,1);
    v1rate_r  = zeros(ln,1);
    Vdot      = zeros(ln,1);
    
end

%% FOI:

phi = 1;%+data.amp*cos((t-32-data.phi)/(365/2*pi));

ddk    = 10^5*sum(mu.*(H+Hv1))/sum(data.Npop);
sd_fun = @(l,b,c,t,d) l + (1-l)*exp(-b*exp(-c*t)*d);%here, t is time since response

if strcmp(data.inp3,'No Closures')||i==1;
    betamod = 1;
elseif any(i==data.imand);
    betamod = min(sd_fun(p2.sdl,p2.sdb,p2.sdc,t-p2.Tres,ddk), sd_fun(p2.sdl,p2.sdb,p2.sdc,14,2));
else
    betamod = sd_fun(p2.sdl,p2.sdb,p2.sdc,t-p2.Tres,ddk);
end

I       = (red*Ina+Ins) + (1-trv1)*(red*Inav1+Insv1) + tm_a*red*(Isa+(1-trv1)*Isav1) + tm_s.*(Iss+(1-trv1)*Issv1);
foi     = phi*beta*betamod*(D*(I./NN));

seedvec = 1e-9*sum(data.Npop)*data.Ev(:,i);%one billionth of the population
seed    = phi*beta*betamod*(D*(seedvec./NN));

%% EQUATIONS:

Sndot    = -Sn.*(foi+seed) - v1rate_sn;

Sdot     = -S.*(foi+seed) + nu.*R - v1rate_sn - v1rate_sw + nuv1.*Sv1;
Shv1dot  = v1rate_sn - hrv1*Shv1 - Shv1.*(foi+seed);
Sv1dot   = hrv1*Shv1 - Sv1.*(1-scv1).*(foi+seed) - nuv1.*Sv1;

Edot     = S.*(foi+seed) + Shv1.*(foi+seed) - (sig1+sig2+sig3+sig4).*E;
Ev1dot   = Sv1.*(1-scv1).*(foi+seed) - (sig1+sig2+sig3+sig4).*Ev1;

Inadot   = sig1.*E - g1.*Ina;               
Insdot   = sig2.*E - (g2+h).*Ins;           
Isadot   = sig3.*E - g1.*Isa;               
Issdot   = sig4.*E - (g2+h).*Iss;           

Inav1dot = sig1.*Ev1 - g1.*Inav1;             
Insv1dot = sig2.*Ev1 - (g2_v1+h_v1).*Insv1;   
Isav1dot = sig3.*Ev1 - g1.*Isav1;             
Issv1dot = sig4.*Ev1 - (g2_v1+h_v1).*Issv1;   

Hdot     = h.*Ins + h.*Iss - (g3+mu).*H;
Hv1dot   = h_v1.*Insv1 + h_v1.*Issv1 - (g3+mu).*Hv1;

Rdot     = g1.*Ina + g1.*Isa + g2.*Ins + g2.*Iss + g3.*H - nu.*R - v1rate_r;
Rv1dot   = g1.*Inav1 + g1*Isav1 + g2_v1.*Insv1 + g2_v1.*Issv1 + g3.*Hv1 + v1rate_sw + v1rate_r;

DEdot    = mu.*H + mu.*Hv1;     

%% OUTPUT:

f= [Sdot;Edot;...
    Inadot;Isadot;Insdot;Issdot;...
    Hdot;Rdot;...
    Shv1dot;Sv1dot;Ev1dot;...
    Inav1dot;Isav1dot;Insv1dot;Issv1dot;...
    Hv1dot;Rv1dot;...
    DEdot;Vdot;Sndot];
f(y<eps) = max(0,f(y<eps)); 

g=h.*(Ins+Iss)+h_v1.*(Insv1+Issv1);%Hin

end

%%

function [value,isterminal,direction] = unmitigated(t,y,data,dis,i,p2)
    
    ln   = length(data.NNs);
    
    S    = y(0*ln+1:1*ln);
    Shv1 = y(8*ln+1:9*ln);
    Sv1  = y(9*ln+1:10*ln);
    Sn   = y(19*ln+1:20*ln);
    
    amp  = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
    ph   = amp.*dis.ph;
    Ts   = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    g2   = (1-ph)./Ts;
    h    = ph./Ts;
    
    %% event 1: end of testing
    %stop testing at first occurence of: Rt<1 if lifted, end of vaccination campaign, 2.5 years after response time
    E1iflag = floor(i/5);
    E1tflag = max(0,data.tvec(end-1)+0.1-t);
    E1vflag = max(0,p2.end-t)*max(0,p2.Tres+2.5*365-t);
    if E1iflag == 0 && E1tflag == 0 && E1vflag ~=0;
        Rt2     = rep_num(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
        E1vflag = max(0,Rt2-1);
    end
    
    value(1)      = E1iflag + E1tflag + E1vflag;
    direction(1)  = -1;
    isterminal(1) = 1;
    
end

function [value,isterminal,direction] = reactive_closures(t,y,data,dis,i,p2)
    
    ln    = length(data.NNs);
    
    S     = y(0*ln+1:1*ln);
    Ins   = y(4*ln+1:5*ln);
    Iss   = y(5*ln+1:6*ln);
    H     = y(6*ln+1:7*ln);
    Shv1  = y(8*ln+1:9*ln);
    Sv1   = y(9*ln+1:10*ln);
    Insv1 = y(13*ln+1:14*ln);
    Issv1 = y(14*ln+1:15*ln);
    Hv1   = y(15*ln+1:16*ln);
    Sn    = y(19*ln+1:20*ln);
    
    occ    = max(1,sum(H+Hv1));
    Hmax   = p2.Hmax;
    amp    = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
    th0    = max(1,1+1.87*((occ-Hmax)/(2*Hmax-Hmax)));
    ph     = amp.*dis.ph;
    pd     = min(th0*dis.pd,1);
    Ts     = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    Th     = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
    g2     = (1-ph)./Ts;
    g3     = (1-pd)./Th;
    h      = ph./Ts;
    mu     = pd./Th;
    h_v1   = dis.h_v1;
    Hdot   = h.*Ins      +h.*Iss        -(g3+mu).*H;
    Hv1dot = h_v1.*Insv1 +h_v1.*Issv1   -(g3+mu).*Hv1;
    occdot = sum(Hdot+Hv1dot);
    r      = occdot/occ;
    Tcap   = t + log(p2.Hmax/occ)/r;
    Tld    = Tcap - 5;%(20 + 5*log(abs(r)));%empirical function, unused for r < 0.025 as below
    
    %% event 1: first measures
    %home-working and distancing imposed at response time
    E1iflag = abs(i-1);
    E1tflag = max(0,data.tvec(end-1)+0.1-t);
    E1vflag = max(0,p2.Tres-t);
        
    value(1)      = E1iflag + E1tflag + E1vflag;
    direction(1)  = -1;
    isterminal(1) = 1;
    
    %% event 2: early lockdown
    %lockdown if less than 4 days before hospital capacity expected to be breached and growth rate is large
    E2iflag = abs((i-2)*(i-4));
    E2tflag = max(0,data.tvec(end-1)+0.1-t);
    E2vflag = max(0,Tld-t) + max(0,0.025-r);  
    
    value(2)      = E2iflag + E2tflag + E2vflag;
    direction(2)  = -1;
    isterminal(2) = 1;
    
    %% event 3: late lockdown
    %lockdown if hospital occupancy greater than 95% of capacity
    E3iflag = abs((i-1)*(i-2)*(i-4));
    E3tflag = max(0,data.tvec(end-1)+0.1-t);
    E3vflag = max(0,0.95*p2.Hmax-occ);  
    
    value(3)      = E3iflag + E3tflag + E3vflag;
    direction(3)  = -1;
    isterminal(3) = 1;
    
    %% event 4: partial reopening
    %partially reopen after 1 week if hospital occupancy less than 25% of capacity
    E4iflag = abs(i-3);
    E4tflag = max(0,data.tvec(end-1)+7-t);
    E4vflag = max(0,occ-0.25*p2.Hmax);  
    
    value(4)      = E4iflag + E4tflag + E4vflag;
    direction(4)  = -1;
    isterminal(4) = 1;
    
    %% event 5: end of closures and testing
    %remove measures at first occurence of: Rt<1 if lifted, end of vaccination campaign and (below 25% occupancy or low non-lockdown growth rate or 90 days since end of rollout), 2.5 years after response time
    E5iflag = floor(i/5);
    E5tflag = max(0,data.tvec(end-1)+0.1-t);
    E5vflag = (max(0,p2.end-t) + max(0,occ-0.25*p2.Hmax)*(max(0,r-0.025) + abs((i-1)*(i-2)*(i-4)))*max(0,p2.end+90-t))*max(0,p2.Tres+2.5*365-t);
    if E5iflag == 0 && E5tflag == 0 && E5vflag ~=0;
        Rt2     = rep_num(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
        E5vflag = max(0,Rt2-1);
    end
    
    value(5)      = E5iflag + E5tflag + E5vflag;
    direction(5)  = -1;
    isterminal(5) = 1;
    
end

function [value,isterminal,direction] = elimination(t,y,data,dis,i,p2)
    
    ln    = length(data.NNs);
    
    S     = y(0*ln+1:1*ln);
    E     = y(1*ln+1:2*ln);
    H     = y(6*ln+1:7*ln);
    Shv1  = y(8*ln+1:9*ln);
    Sv1   = y(9*ln+1:10*ln);
    Ev1   = y(10*ln+1:11*ln);
    Hv1   = y(15*ln+1:16*ln);
    Sn    = y(19*ln+1:20*ln);
    
    occ   = max(1,sum(H+Hv1));
    amp   = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
    ph    = amp.*dis.ph;
    Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
    g2    = (1-ph)./Ts;
    h     = ph./Ts;
    if t<p2.t_tit;   
        asc_a = 0;
        asc_s = 0;
        tm_a  = 1;
        tm_s  = 1;
        
    elseif t<p2.end && i~=5;
        trate  = p2.trate;
        incid  = max(0,10^5*((dis.siga+dis.sigs)*sum(E+Ev1))/sum(data.Npop));
        asc_s  = 1/(1+exp(0.7623+1.605*log10(incid)-1.416*log10(trate)));
        propCT = 1/(1+exp(2.159+1.697*log10(incid)));
        asc_a  = propCT*asc_s + (1-propCT)*0;
        
        asc_a  = min(asc_a,(trate/incid)*(0*(1-propCT) + (1-dis.ps)*propCT));
        asc_s  = min(asc_s,(trate/incid)*(1*(1-propCT) + dis.ps*propCT));
        asc_a  = max(trate/10^5*(0*(1-propCT) + (1-dis.ps)*propCT),asc_a);
        asc_s  = max(trate/10^5*(1*(1-propCT) + dis.ps*propCT),asc_s);
        
        onsPCR_s = 11.3224-2.6260*log10(trate);
        onsPCR_c = onsPCR_s - 5.6304;
        Teff_c   = max(0,dis.Tinc+onsPCR_c-dis.Tlat);
        Teff_s   = max(0,dis.Tinc+onsPCR_s-dis.Tlat);
        mult_ac  = min(Teff_c,dis.Tay)./dis.Tay;
        mult_sc  = min(Teff_c,Ts)./Ts;
        mult_ss  = min(Teff_s,Ts)./Ts;
        
        tm_a = mult_ac;
        tm_s = mult_sc*propCT + mult_ss*(1-propCT);
    
    else       
        asc_a = 0;
        asc_s = 0;
        tm_a  = 1;
        tm_s  = 1;
            
    end
    sig1  = dis.siga*(1-asc_a);
    sig2  = dis.sigs*(1-asc_s);
    sig3  = dis.siga*asc_a;
    sig4  = dis.sigs*asc_s;
    
    %% event 1: first lockdown
    %lockdown at first occurence of: response time, 95% of hospital capacity
    E1iflag = abs(i-1);
    E1tflag = max(0,data.tvec(end-1)+0.1-t);
    E1vflag = max(0,p2.Tres-t)*max(0,0.95*p2.Hmax-occ);
        
    value(1)      = E1iflag + E1tflag + E1vflag;
    direction(1)  = -1;
    isterminal(1) = 1;
    
    %% event 2: domestic reopening
    %reopen domestic economy after 1 week if Rt<1  
    E2iflag = abs(i-2);
    E2tflag = max(0,data.tvec(end-1)+7-t);
    E2vflag = 1;
    if E2iflag == 0 && E2tflag == 0;
        Rt1     = rep_num(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,3),1,sig1,sig2,sig3,sig4,tm_a,tm_s);
        E2vflag = max(0,Rt1-1);
    end

    value(2)      = E2iflag + E2tflag + E2vflag;
    direction(2)  = -1;
    isterminal(2) = 1;
    
    %% event 3: relockdown
    %lockdown again after 1 week if Rt>1.2
    E3iflag = abs(i-3);
    E3tflag = max(0,data.tvec(end-1)+7-t);
    E3vflag = 1;
    if E3iflag == 0 && E3tflag == 0;
        Rt1     = rep_num(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,3),1,sig1,sig2,sig3,sig4,tm_a,tm_s);
        E3vflag = max(0,1.2-Rt1);
    end

    value(3)      = E3iflag + E3tflag + E3vflag;
    direction(3)  = -1;
    isterminal(3) = 1;
    
    %% event 4: end of closures and testing
    %remove measures at first occurence of: Rt<1 if lifted, end of vaccination campaign, 2.5 years after response time
    E4iflag = floor(i/5);
    E4tflag = max(0,data.tvec(end-1)+0.1-t);
    E4vflag = max(0,p2.end-t)*max(0,p2.Tres+2.5*365-t);
    if E4iflag == 0 && E4tflag == 0 && E4vflag ~=0;
        Rt2     = rep_num(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
        E4vflag = max(0,Rt2-1);
    end
    
    value(4)      = E4iflag + E4tflag + E4vflag;
    direction(4)  = -1;
    isterminal(4) = 1;
    
end

function Rt = rep_num(dis,h,g2,S,Shv1,Sv1,N,D,betamod,sig1,sig2,sig3,sig4,tm_a,tm_s)
    
    N(N==0) = 1;
    ln      = length(N);
    F       = zeros(10*ln,10*ln);
    
    FOIu = dis.beta.*betamod.*D./repmat(N',ln,1).*repmat(S+Shv1,1,ln);
    FOIv = dis.beta.*betamod.*(1-dis.scv1).*D./repmat(N',ln,1).*repmat(Sv1,1,ln);
    
    F(1:ln,     2*ln+1:end) = [dis.red.*FOIu,  FOIu,  tm_a.*dis.red*FOIu,  tm_s.*FOIu, ...  
                               dis.red.*(1-dis.trv1).*FOIu,  (1-dis.trv1).*FOIu,  tm_a.*dis.red.*(1-dis.trv1).*FOIu,  tm_s.*(1-dis.trv1).*FOIu];
    F(ln+1:2*ln,2*ln+1:end) = [dis.red.*FOIv,  FOIv,  tm_a.*dis.red.*FOIv,  tm_s.*FOIv, ...  
                               dis.red.*(1-dis.trv1).*FOIv,  (1-dis.trv1).*FOIv,  tm_a.*dis.red.*(1-dis.trv1).*FOIv,  tm_s.*(1-dis.trv1).*FOIv];

    onesn = ones(ln,1);
    vvec  = [(sig1+sig2+sig3+sig4).*onesn;
             (sig1+sig2+sig3+sig4).*onesn;
              dis.g1.*onesn;
             (g2+h).*onesn;
              dis.g1.*onesn;
             (g2+h).*onesn;
              dis.g1.*onesn;
             (dis.g2_v1+dis.h_v1).*onesn;
              dis.g1.*onesn;
             (dis.g2_v1+dis.h_v1).*onesn];
    V     = diag(vvec);
    
    V(2*ln+1:3*ln,1:ln)       = diag(-sig1.*onesn);
    V(3*ln+1:4*ln,1:ln)       = diag(-sig2.*onesn);
    V(4*ln+1:5*ln,1:ln)       = diag(-sig3.*onesn);
    V(5*ln+1:6*ln,1:ln)       = diag(-sig4.*onesn);
    V(6*ln+1:7*ln,ln+1:2*ln)  = diag(-sig1.*onesn);
    V(7*ln+1:8*ln,ln+1:2*ln)  = diag(-sig2.*onesn);
    V(8*ln+1:9*ln,ln+1:2*ln)  = diag(-sig3.*onesn);
    V(9*ln+1:10*ln,ln+1:2*ln) = diag(-sig4.*onesn);
    
    NGM = F/V;
    Rt  = eigs(NGM,1);
    
end
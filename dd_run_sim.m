function [data,f,g] = dd_run_sim(data,dis,p2);

ln    = length(data.NNs);
lx    = length(data.obj);
adInd = 3;

t0   = data.tvec(1);
y0   = dis.y0;
i    = 1;
tend = data.tvec(end);

output = struct;

while i < 6; 
    D = data.Dvec(:,:,i);
    
    fun                = @(t,y)ODEs(data,D,i,t,dis,y,p2);
    options            = odeset('Events', @(t,y) data.ev_fn(t,y,data,dis,i,p2));%,'MaxStep',0.1);   
    [tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);
     output            = dd_run_postprocessing(data,dis,i,p2,output,tout,yout);

    if ~isempty(ie);
        inext = data.inext(ie(end));
    else
        inext = 6;
    end

    if inext < 6;
        t0        = tout(end);
        y0        = yout(end,:)';
        i         = inext;
        data.tvec = [data.tvec(1:end-1),tout(end),tend];
        p2.Tres   = min(p2.Tres,tout(end));
    else
        i              = inext;
        data.tvec(end) = tout(end);
        f              = output.f;
        g              = output.g;
    end    
end

end

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

occ   = sum(H+Hv1);
Hmax  = p2.Hmax;
th    = p2.th;

%% TIME-DEPENDENT DISEASE PARAMETERS:

%Amplitudes
amp = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
th0 = max(1,1+th*((occ-Hmax)/(2*Hmax-Hmax)));

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
asca    = p2.asca;
ascb    = p2.ascb;
ascc    = p2.ascc;
pcta    = p2.pcta;
pctb    = p2.pctb;
opsa    = p2.opsa;
opsb    = p2.opsb;
opc     = p2.opc;
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
    asc_s  = 1/(1+exp(asca+ascb*log10(incid)+ascc*log10(trate)));
    propCT = 1/(1+exp(pcta+pctb*log10(incid)));
    asc_a  = propCT*asc_s + (1-propCT)*0;
    
    asc_a = min(asc_a,(trate/incid)*(0*(1-propCT) + (1-ps)*propCT));
    asc_s = min(asc_s,(trate/incid)*(1*(1-propCT) + ps*propCT));
    asc_a = max(trate/10^5*(0*(1-propCT) + (1-ps)*propCT),asc_a);
    asc_s = max(trate/10^5*(1*(1-propCT) + ps*propCT),asc_s);
    
    onsPCR_s = opsa+opsb*log10(trate);
    onsPCR_c = onsPCR_s+opc;
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

%the vaccine rollout is not perfectly efficient, only susceptibles benefit from vaccination
nonVax = S + E + Ina + Ins + Isa + Iss + H + R + DE;

if t>=pend
    vrate_s = zeros(ln,1);
    Vdot    = zeros(ln,1);    
elseif t>=startp4
    vrate_s = aratep4.*S./nonVax;
    Vdot    = aratep4;
elseif t>=startp3
    vrate_s = aratep3.*S./nonVax;
    Vdot    = aratep3;
elseif t>=startp2
    vrate_s = aratep2.*S./nonVax;
    Vdot    = aratep2;
elseif t>=startp1
    vrate_s = aratep1.*S./nonVax;
    Vdot    = aratep1;
else
    vrate_s = zeros(ln,1);
    Vdot    = zeros(ln,1);
end

%% FOI:

phi = 1;%+data.amp*cos((t-32-data.phi)/(365/2*pi));

ddk    = max(0,10^5*sum(mu.*(H+Hv1))/sum(data.Npop));
sd_fun = @(a,b,c,t,d) 1/(1 + exp(a + b*log10(d) - c*t));%here, t is time since response

if strcmp(data.inp3,'No Closures')||i==1;
    betamod = 1;
elseif any(i==data.imand);
    betamod = min(sd_fun(p2.sda,p2.sdb,p2.sdc,t-p2.Tres,ddk), sd_fun(p2.sda,p2.sdb,p2.sdc,14,2));
else
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,t-p2.Tres,ddk);
end

I       = (red*Ina+Ins) + (1-trv1)*(red*Inav1+Insv1) + tm_a*red*(Isa+(1-trv1)*Isav1) + tm_s.*(Iss+(1-trv1)*Issv1);
foi     = phi*beta*betamod*(D*(I./NN));

seedvec = 1e-9*sum(data.Npop)*dis.Ev*data.xconf(i,data.IntlInd);%one billionth of the population
seed    = phi*beta*betamod*(D*(seedvec./NN));

%% EQUATIONS:

Sndot    = -Sn.*(foi+seed) - vrate_s;

Sdot     = -S.*(foi+seed) + nu.*R - vrate_s + nuv1.*Sv1;%- v1rate_sw
Shv1dot  = vrate_s - hrv1*Shv1 - Shv1.*(foi+seed);
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

Rdot     = g1.*Ina + g2.*Ins + g1.*Isa + g2.*Iss + g3.*H - nu.*R;% - v1rate_r;
Rv1dot   = g1.*Inav1 + g2_v1.*Insv1 + g1*Isav1 + g2_v1.*Issv1 + g3.*Hv1;% + v1rate_r + v1rate_sw

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
function [p2,f,g] = dd_run_sim(data,dis,p2)

t0      = data.tvec(1);
y0      = dis.y0;
i       = 1;
tend    = data.tvec(end);
p2.tvec = data.tvec;
output  = struct;

while i < 6; 
    D = data.Dvec(:,:,i);
    
    fun                = @(t,y) ODEs(t,y,data,dis,i,D,p2);
    options            = odeset('Events', @(t,y) data.ev_fn(t,y,data,dis,i,D,p2));%,'MaxStep',0.1);   
    [tout,yout,~,~,ie] = ode45(fun,[t0 tend],y0,options);
    
    if ~isempty(ie);
        inext = data.inext(ie(end));
    else
        inext = 6;
    end
    output = dd_run_postprocessing(data,dis,i,p2,output,tout,yout,inext);
    
    if inext < 6;
        t0        = tout(end);
        y0        = yout(end,:)';
        i         = inext;
        p2.tvec   = [p2.tvec(1:end-1),t0,tend];
        p2.Tres   = min(p2.Tres,tout(end));
        p2.tmand  = t0-p2.Tres;
        p2.dmand  = output.dmand;
    else
        i              = inext;
        data.tvec(end) = tout(end);
        f              = output.f;
        g              = output.g;
    end    
end

end

function dydt = ODEs(t,y,data,dis,i,D,p2)

ln = length(data.NNs);

%% INITIAL CONDITION

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

%% VACCINATION

%the vaccine rollout is not perfectly efficient, only susceptibles benefit from vaccination
nonVax = S + E + Ina + Ins + Isa + Iss + H + R + DE;

if t >= p2.end;
    vrate_s = zeros(ln,1);
    Vdot    = zeros(ln,1);    
elseif t >= p2.startp4;
    vrate_s = p2.aratep4.*S./nonVax;
    Vdot    = p2.aratep4;
elseif t >= p2.startp3;
    vrate_s = p2.aratep3.*S./nonVax;
    Vdot    = p2.aratep3;
elseif t >= p2.startp2;
    vrate_s = p2.aratep2.*S./nonVax;
    Vdot    = p2.aratep2;
elseif t >= p2.startp1;
    vrate_s = p2.aratep1.*S./nonVax;
    Vdot    = p2.aratep1;
else
    vrate_s = zeros(ln,1);
    Vdot    = zeros(ln,1);
end

amp = min((Sn + (S-Sn).*(1-dis.heff))./S, 1);%vaccine waning reduces ph
ph  = amp.*dis.ph;
Ts  = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g2  = (1-ph)./Ts;
h   = ph./Ts;

%% HEALTHCARE

occ = sum(H+Hv1);
th0 = max(1, 1+p2.th*((occ-p2.Hmax)/p2.Hmax));%overcapacity hospitals increases pd
pd  = min(th0*dis.pd,1);
Th  = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
g3  = (1-pd)./Th;
mu  = pd./Th;

%% CASE ISOLATION AND TRACING

cit_sw  = double((t >= p2.t_tit) & (i~=5));
prev_sw = double(sum(E + Ev1 + Ina + Isa + Ins + Iss + Inav1 + Isav1 + Insv1 + Issv1 + H + Hv1) < 10^-7*sum(data.Npop));

Tss_eff = max(0,p2.iisym-dis.Tlat);
Tsa_eff = max(0,p2.iitra-dis.Tlat);
tm_a    = 1*(1-cit_sw*prev_sw) + (min(Tsa_eff,dis.Tay)./dis.Tay)*(cit_sw*prev_sw);
tm_s    = 1*(1-cit_sw) + (min(Tss_eff,Ts)./Ts)*(cit_sw*(1-prev_sw)) + (min(Tsa_eff,Ts)./Ts)*(cit_sw*prev_sw);

% DISTANCING
ddk    = max(0,10^5*sum(mu.*(H+Hv1))/sum(data.Npop));
sd_fun = @(a,b,c,t,d) 1/(1 + exp(a + b*log10(d) - c*t));%here, t is time since response time
if i==1;%strcmp(data.inp3,'No Closures')||
    betamod = 1;
elseif any(i==data.imand);
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,p2.tmand,max(p2.dmand,ddk));
else
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,t-p2.Tres,ddk);
end
I         = (dis.red*Ina+Ins) + (1-dis.trv1)*(dis.red*Inav1+Insv1) + tm_a*dis.red*(Isa+(1-dis.trv1)*Isav1) + tm_s.*(Iss+(1-dis.trv1)*Issv1);
NN        = data.NNs;
NN(NN==0) = 1;
foi       = dis.beta*betamod*(D*(I./NN));
seedvec   = 1e-8*sum(data.Npop)*dis.Ev*data.xconf(i,data.IntlInd);%one billionth of the population
seed      = dis.beta*betamod*(D*(seedvec./NN));

inc   = sum(S.*(foi+seed) + Shv1.*(foi+seed) + Sv1.*(1-dis.scv1).*(foi+seed));
sinc  = dis.ps*inc;
hadm  = sum(h.*Ins + h.*Iss + dis.h_v1.*Insv1 + dis.h_v1.*Issv1);
asc_a = 0*(1-cit_sw*prev_sw) + min(p2.masc, max(0, p2.trate - hadm)/inc)*(cit_sw*prev_sw);
asc_s = 0*(1-cit_sw) + min(p2.masc, max(0, p2.trate - hadm)/sinc)*(cit_sw*(1-prev_sw)) + min(p2.masc, max(0, p2.trate - hadm)/inc)*(cit_sw*prev_sw);
sig1  = dis.siga*(1-asc_a);
sig2  = dis.sigs*(1-asc_s);
sig3  = dis.siga*asc_a;
sig4  = dis.sigs*asc_s; 

%% EQUATIONS

Sdot     = -S.*(foi+seed) + dis.nu.*R - vrate_s + dis.nuv1.*Sv1;%- v1rate_sw
Shv1dot  = vrate_s - dis.hrv1*Shv1 - Shv1.*(foi+seed);
Sv1dot   = dis.hrv1*Shv1 - Sv1.*(1-dis.scv1).*(foi+seed) - dis.nuv1.*Sv1;

Edot     = S.*(foi+seed) + Shv1.*(foi+seed) - (sig1+sig2+sig3+sig4).*E;
Ev1dot   = Sv1.*(1-dis.scv1).*(foi+seed) - (sig1+sig2+sig3+sig4).*Ev1;

Inadot   = sig1.*E - dis.g1.*Ina;               
Insdot   = sig2.*E - (g2+h).*Ins;           
Isadot   = sig3.*E - dis.g1.*Isa;               
Issdot   = sig4.*E - (g2+h).*Iss;           
Inav1dot = sig1.*Ev1 - dis.g1.*Inav1;             
Insv1dot = sig2.*Ev1 - (dis.g2_v1+dis.h_v1).*Insv1;   
Isav1dot = sig3.*Ev1 - dis.g1.*Isav1;             
Issv1dot = sig4.*Ev1 - (dis.g2_v1+dis.h_v1).*Issv1;   

Hdot     = h.*Ins + h.*Iss - (g3+mu).*H;
Hv1dot   = dis.h_v1.*Insv1 + dis.h_v1.*Issv1 - (g3+mu).*Hv1;

Rdot     = dis.g1.*Ina + g2.*Ins + dis.g1.*Isa + g2.*Iss + g3.*H - dis.nu.*R;% - v1rate_r;
Rv1dot   = dis.g1.*Inav1 + dis.g2_v1.*Insv1 + dis.g1*Isav1 + dis.g2_v1.*Issv1 + g3.*Hv1;% + v1rate_r + v1rate_sw

DEdot    = mu.*H + mu.*Hv1;     

Sndot    = -Sn.*(foi+seed) - vrate_s.*(Sn./S);

%% OUTPUT

dydt        = [Sdot;Edot;Inadot;Isadot;Insdot;Issdot;Hdot;Rdot;Shv1dot;Sv1dot;Ev1dot;Inav1dot;Isav1dot;Insv1dot;Issv1dot;Hv1dot;Rv1dot;DEdot;Vdot;Sndot];
dydt(y<eps) = max(0,dydt(y<eps)); 

end
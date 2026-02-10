function [value,isterminal,direction] = dd_set_rules_elimination(t,y,data,dis,i,D,p2)
    
ln        = length(data.NNs);

S         = y(0*ln+1:1*ln);
E         = y(1*ln+1:2*ln);
Ina       = y(2*ln+1:3*ln);
Isa       = y(3*ln+1:4*ln);
Ins       = y(4*ln+1:5*ln);
Iss       = y(5*ln+1:6*ln);
H         = y(6*ln+1:7*ln);
R         = y(7*ln+1:8*ln);
Shv1      = y(8*ln+1:9*ln);
Sv1       = y(9*ln+1:10*ln);
Ev1       = y(10*ln+1:11*ln);
Inav1     = y(11*ln+1:12*ln);
Isav1     = y(12*ln+1:13*ln);
Insv1     = y(13*ln+1:14*ln);
Issv1     = y(14*ln+1:15*ln);
Hv1       = y(15*ln+1:16*ln);
Rv1       = y(16*ln+1:17*ln);
DE        = y(17*ln+1:18*ln);
Sn        = y(19*ln+1:20*ln);

amp       = min((Sn + (S-Sn).*(1-dis.heff))./S, 1);
ph        = amp.*dis.ph;
Ts        = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g2        = (1-ph)./Ts;
h         = ph./Ts;

occ       = sum(H+Hv1);
th0       = max(1, 1+p2.th*((occ-p2.Hmax)/p2.Hmax));%overcapacity hospitals increases pd
pd        = min(th0*dis.pd,1);
Th        = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
mu        = pd./Th;

cit_sw    = double((t >= p2.t_tit) & (i~=5));
prev_sw   = double(sum(E + Ev1 + Ina + Isa + Ins + Iss + Inav1 + Isav1 + Insv1 + Issv1 + H + Hv1) < 10^-7*sum(data.Npop));

Tss_eff   = max(0,p2.iisym-dis.Tlat);
Tsa_eff   = max(0,p2.iitra-dis.Tlat);
tm_a      = 1*(1-cit_sw*prev_sw) + (min(Tsa_eff,dis.Tay)./dis.Tay)*(cit_sw*prev_sw);
tm_s      = 1*(1-cit_sw) + (min(Tss_eff,Ts)./Ts)*(cit_sw*(1-prev_sw)) + (min(Tsa_eff,Ts)./Ts)*(cit_sw*prev_sw);

% DISTANCING
ddk       = max(0,10^5*sum(mu.*(H+Hv1))/sum(data.Npop));
sd_fun    = @(a,b,c,t,d) 1/(1 + exp(a + b*log10(d) - c*t));%here, t is time since response time
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

inc       = sum(S.*(foi+seed) + Shv1.*(foi+seed) + Sv1.*(1-dis.scv1).*(foi+seed));
sinc      = dis.ps*inc;
hadm      = sum(h.*Ins + h.*Iss + dis.h_v1.*Insv1 + dis.h_v1.*Issv1);
asc_a     = 0*(1-cit_sw*prev_sw) + min(p2.masc, max(0, p2.trate - hadm)/inc)*(cit_sw*prev_sw);
asc_s     = 0*(1-cit_sw) + min(p2.masc, max(0, p2.trate - hadm)/sinc)*(cit_sw*(1-prev_sw)) + min(p2.masc, max(0, p2.trate - hadm)/inc)*(cit_sw*prev_sw);
sig1      = dis.siga*(1-asc_a);
sig2      = dis.sigs*(1-asc_s);
sig3      = dis.siga*asc_a;
sig4      = dis.sigs*asc_s; 

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
E2tflag = max(0,data.tvec(end-1)+14-t);
E2vflag = 1;
if E2iflag == 0 && E2tflag == 0;
    [Rt1,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,3),1,sig1,sig2,sig3,sig4,tm_a,tm_s);
    E2vflag = max(0,Rt1-1);
end

value(2)      = E2iflag + E2tflag + E2vflag;
direction(2)  = -1;
isterminal(2) = 1;

%% event 3: relockdown
%lockdown again after 1 week if Rt>1.2
E3iflag = abs(i-3);
E3tflag = max(0,data.tvec(end-1)+14-t);
E3vflag = 1;
if E3iflag == 0 && E3tflag == 0;
    [Rt1,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,3),1,sig1,sig2,sig3,sig4,tm_a,tm_s);
    E3vflag = max(0,1.2-Rt1);
end

value(3)      = E3iflag + E3tflag + E3vflag;
direction(3)  = -1;
isterminal(3) = 1;

%% event 4: end of closures and testing
%remove measures at first occurence of: Rt<1 if lifted, end of vaccination campaign, 2.5 years after response time
E4iflag = abs((i-2)*(i-3));
E4tflag = max(0,data.tvec(end-1)+0.1-t);
E4vflag = max(0,p2.end-t)*max(0,p2.Tres+2.5*365-t);
if E4iflag == 0 && E4tflag == 0 && E4vflag ~=0;
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E4vflag = max(0,Rt2-1);
end

value(4)      = E4iflag + E4tflag + E4vflag;
direction(4)  = -1;
isterminal(4) = 1;

%% event 5: end of simulation
%stop simulation when all measures have been removed, the DFE has been reached and Rt<=1
E5iflag = abs(i-5);
E5tflag = max(0,data.tvec(end-1)+0.1-t);
E5vflag = max(0,0.99*sum(data.NNs)-sum(S+R+Rv1+DE)) + max(0,occ-0.000001*sum(data.NNs));
if E5iflag + E5tflag + E5vflag == 0;
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E5vflag = max(0,Rt2-1);
end

value(5)      = E5iflag + E5tflag + E5vflag;
direction(5)  = -1;
isterminal(5) = 1;
    
end
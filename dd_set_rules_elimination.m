function [value,isterminal,direction] = dd_set_rules_elimination(t,y,data,dis,i,p2)
    
ln    = length(data.NNs);

S     = y(0*ln+1:1*ln);
E     = y(1*ln+1:2*ln);
H     = y(6*ln+1:7*ln);
R     = y(7*ln+1:8*ln);
Shv1  = y(8*ln+1:9*ln);
Sv1   = y(9*ln+1:10*ln);
Ev1   = y(10*ln+1:11*ln);
Hv1   = y(15*ln+1:16*ln);
Sn    = y(19*ln+1:20*ln);
Rv1   = y(16*ln+1:17*ln);
DE    = y(17*ln+1:18*ln);

amp   = min((Sn + (S-Sn).*(1-dis.heff))./S, 1);
ph    = amp.*dis.ph;
Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g2    = (1-ph)./Ts;
h     = ph./Ts;
occ   = sum(H+Hv1);
if t>=p2.t_tit && i~=5;
    incid  = max(0,10^5*((dis.siga+dis.sigs)*sum(E+Ev1))/sum(data.Npop));
    asc_s  = 1/(1+exp(p2.asca + p2.ascb*log10(incid) + p2.ascc*log10(p2.trate)));
    propCT = 1/(1+exp(p2.pcta + p2.pctb*log10(incid)));
    asc_a  = propCT*asc_s + (1-propCT)*0;
    
    asc_a = min(asc_a,(p2.trate/incid)*(0*(1-propCT) + (1-dis.ps)*propCT));
    asc_s = min(asc_s,(p2.trate/incid)*(1*(1-propCT) + dis.ps*propCT));
    asc_a = max(p2.trate/10^5*(0*(1-propCT) + (1-dis.ps)*propCT),asc_a);
    asc_s = max(p2.trate/10^5*(1*(1-propCT) + dis.ps*propCT),asc_s);
    
    onsPCR_s = p2.opsa + p2.opsb*log10(p2.trate);
    onsPCR_c = onsPCR_s + p2.opc;
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

sig1 = dis.siga*(1-asc_a);
sig2 = dis.sigs*(1-asc_s);
sig3 = dis.siga*asc_a;
sig4 = dis.sigs*asc_s;

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
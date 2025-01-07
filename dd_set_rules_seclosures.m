function [value,isterminal,direction] = dd_set_rules_seclosures(t,y,data,dis,i,p2)
    
ln     = length(data.NNs);

S      = y(0*ln+1:1*ln);
Ins    = y(4*ln+1:5*ln);
Iss    = y(5*ln+1:6*ln);
H      = y(6*ln+1:7*ln);
R      = y(7*ln+1:8*ln);
Shv1   = y(8*ln+1:9*ln);
Sv1    = y(9*ln+1:10*ln);
Insv1  = y(13*ln+1:14*ln);
Issv1  = y(14*ln+1:15*ln);
Hv1    = y(15*ln+1:16*ln);
Sn     = y(19*ln+1:20*ln);
Rv1    = y(16*ln+1:17*ln);
DE     = y(17*ln+1:18*ln);

occ    = max(0.01,sum(H+Hv1));%capped due to division by occ below
Hmax   = p2.Hmax;
th     = p2.th;
amp    = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
th0    = max(1,1+th*((occ-Hmax)/(2*Hmax-Hmax)));
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
Tgen   = dis.Tlat + dis.Tsh;
Tld    = Tcap - Tgen/2;%empirical function, unused for r < 0.025 as below

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
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E5vflag = max(0,Rt2-1);
end

value(5)      = E5iflag + E5tflag + E5vflag;
direction(5)  = -1;
isterminal(5) = 1;

%% event 6: end of simulation
%stop simulation when all measures have been removed, the DFE has been reached and Rt<=1
E6iflag = abs(i-5);
E6tflag = max(0,data.tvec(end-1)+0.1-t);
E6vflag = max(0,0.99*sum(data.NNs)-sum(S+R+Rv1+DE)) + max(0,occ-0.000001*sum(data.NNs));
if E6iflag + E6tflag + E6vflag == 0;
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E6vflag = max(0,Rt2-1);
end

value(6)      = E6iflag + E6tflag + E6vflag;
direction(6)  = -1;
isterminal(6) = 1;
    
end
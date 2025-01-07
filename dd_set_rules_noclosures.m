function [value,isterminal,direction] = dd_set_rules_noclosures(t,y,data,dis,i,p2)
    
ln   = length(data.NNs);

S    = y(0*ln+1:1*ln);
H    = y(6*ln+1:7*ln);
R    = y(7*ln+1:8*ln);
Shv1 = y(8*ln+1:9*ln);
Sv1  = y(9*ln+1:10*ln);
Hv1  = y(15*ln+1:16*ln);
Rv1  = y(16*ln+1:17*ln);
DE   = y(17*ln+1:18*ln);
Sn   = y(19*ln+1:20*ln);

amp  = min((Sn + (S-Sn).*(1-dis.heff))./S, 1);
ph   = amp.*dis.ph;
Ts   = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g2   = (1-ph)./Ts;
h    = ph./Ts;
occ  = sum(H+Hv1);

%% event 1: end of testing
%stop testing at first occurence of: Rt<1 if lifted, end of vaccination campaign, 2.5 years after response time
E1iflag = floor(i/5);
E1tflag = max(0,data.tvec(end-1)+0.1-t);
E1vflag = max(0,p2.end-t)*max(0,p2.Tres+2.5*365-t);
if E1iflag == 0 && E1tflag == 0 && E1vflag ~=0;
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E1vflag = max(0,Rt2-1);
end

value(1)      = E1iflag + E1tflag + E1vflag;
direction(1)  = -1;
isterminal(1) = 1;

%% event 2: end of simulation
%stop simulation when all measures have been removed, the DFE has been reached and Rt<=1
E2iflag = abs(i-5);
E2tflag = max(0,data.tvec(end-1)+0.1-t);
E2vflag = max(0,0.99*sum(data.NNs)-sum(S+R+Rv1+DE)) + max(0,occ-0.000001*sum(data.NNs));
if E2iflag + E2tflag + E2vflag == 0;
    [Rt2,~] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,data.NNs,data.Dvec(:,:,5),1,dis.siga,dis.sigs,0,0,1,1);
    E2vflag = max(0,Rt2-1);
end

value(2)      = E2iflag + E2tflag + E2vflag;
direction(2)  = -1;
isterminal(2) = 1;
    
end
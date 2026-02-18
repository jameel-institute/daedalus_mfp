function rt = calc_r0rt(data,dis,p2,k)

if k == 1;
    ik = 3;
elseif k == 2;
    ik = 4;
elseif k == 3;
    ik = 3;
elseif k == 4;
    ik = 4;
elseif k == 5;
    ik = 3;
end
D  = data.Dvec(:,:,ik);
ln = length(D);

Ts      = dis.Ts;
Tss_eff = max(0,p2.iisym-dis.Tlat);
Tsa_eff = max(0,p2.iitra-dis.Tlat);
if k < 5;
    tm_a  = 1;
    tm_s  = min(Tss_eff,Ts)./Ts;
    asc_a = 0;
    asc_s = min(p2.masc,max(0,p2.trate)/eps);
else
    tm_a  = min(Tsa_eff,dis.Tay)./dis.Tay;
    tm_s  = min(Tsa_eff,Ts)./Ts;
    asc_a = min(p2.masc,max(0,p2.trate)/eps);
    asc_s = asc_a;

end
sig1    = dis.siga*(1-asc_a);
sig2    = dis.sigs*(1-asc_s);
sig3    = dis.siga*asc_a;
sig4    = dis.sigs*asc_s;

[Rt,~]  = dd_calc_Rt(dis,dis.h,dis.g2,data.NNs,zeros(ln,1),zeros(ln,1),data.NNs,D,1,sig1,sig2,sig3,sig4,tm_a,tm_s);
rt      = Rt;                    

end
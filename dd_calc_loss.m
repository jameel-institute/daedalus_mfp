function [cost,ccost_t,c] = dd_calc_loss(data,dis,g)

t  = g(:,1);
lx = length(data.obj);
ln = lx+4;

%% VLYL

deaths    = g(end,(1+2*lx+6*ln+1):(1+2*lx+7*ln));
cost(1,:) = deaths;

lyl       = deaths.*dis.lg;
cost(2,:) = lyl;

vlyl      = lyl*data.vly;
cost(3,:) = vlyl;

deaths          = g(:,(1+2*lx+6*ln+1):(1+2*lx+7*ln));
ccost_t(:,1:ln) = deaths.*dis.lg.*data.vly;

%% VSYL

Stu              = lx+2;
students         = data.NNs(Stu);
cost(4,lx+[1,2]) = students;%student numbers not used, but storing supply costs in lx+1 and demand costs in lx+2

%student supply: learning losses due to illness, when schools are open
%assume that the symptomatic period is the same duration as the infectious period (though start/end points may be different)
%assume that the isolation period is the latent period plus maximum infectious period
%not scaled by ability for sfh because the valuation used here accounts for this
%assume isolating symptomatics who are hospitalised isolate for same duration as infectious period
ph           = g(:,1+2*lx+8*ln+1+2*ln+Stu);
phv          = (1-dis.hv1)*dis.ph(Stu);
Ts           = g(:,1+2*lx+8*ln+1+Stu);
Tiso         = dis.Tlat + max([dis.Tay,dis.Tsr,dis.Tsh]);
isoasyu      = (g(:,1+2*lx+0*ln+Stu)./dis.Tay).*Tiso;%.*(1-(1/3));
isoasyv      = (g(:,1+2*lx+1*ln+Stu)./dis.Tay).*Tiso;%.*(1-(1/3));
isosymu      = (g(:,1+2*lx+2*ln+Stu)./Ts).*(Tiso.*(1-ph) + Ts.*ph);              
isosymv      = (g(:,1+2*lx+3*ln+Stu)./dis.Ts_v1(Stu)).*(Tiso.*(1-phv) + dis.Ts_v1(Stu).*phv);
nissym       = g(:,1+2*lx+4*ln+Stu);
hospts       = g(:,1+2*lx+5*ln+Stu);
deaths       = g(:,1+2*lx+6*ln+Stu);
abs_ill      = isoasyu + isoasyv + isosymu + isosymv + nissym + hospts + deaths;%number of students absent
abs_ill_open = abs_ill.*g(:,1+data.EdInd);%number of students absent in open schools
abs_ilop_int = trapz(t,abs_ill_open)/365;%number of student-years of absence
cost(5,lx+1) = abs_ilop_int;
vsyl_sts     = abs_ilop_int*data.vsy;
cost(6,lx+1) = vsyl_sts;

%student demand: learning losses due to school closures
abs_clos     = students.*(1-g(:,1+data.EdInd));%.*(1-(1/3));%number of students absent
abs_clos_int = trapz(t,abs_clos)/365;%(diff(t)'*presl)/365;%number of student-years of absence
cost(5,lx+2) = abs_clos_int;
vsyl_std     = abs_clos_int*data.vsy;
cost(6,lx+2) = vsyl_std;

ccost_t(:,ln+lx+1) = cumtrapz(t,abs_ill_open,1)./365.*data.vsy;
ccost_t(:,ln+lx+2) = cumtrapz(t,abs_clos,1)./365.*data.vsy;

%% GDPL

notEd = [1:(data.EdInd-1),(data.EdInd+1):lx];

%labour supply: GDP loss due to illness, when businesses are open
%assume that the symptomatic period is the same duration as the infectious period (though start/end points may be different)
%assume that the isolation period is the latent period plus maximum infectious period
%assume isolating asymptomatics wfh if capable
%assume isolating symptomatics who are hospitalised isolate for same duration as infectious period
ph            = g(:,1+2*lx+8*ln+1+2*ln+notEd);
phv           = (1-dis.hv1).*dis.ph(notEd)';
Ts            = g(:,1+2*lx+8*ln+1+notEd);
Tiso          = dis.Tlat + max([dis.Tay,dis.Tsr,dis.Tsh]);
hw            = g(:,1+1*lx+notEd);
isoasyu       = (g(:,1+2*lx+0*ln+notEd)./dis.Tay).*Tiso.*(1-hw);
isoasyv       = (g(:,1+2*lx+1*ln+notEd)./dis.Tay).*Tiso.*(1-hw);
isosymu       = (g(:,1+2*lx+2*ln+notEd)./Ts).*(Tiso.*(1-ph) + Ts.*ph);              
isosymv       = (g(:,1+2*lx+3*ln+notEd)./dis.Ts_v1(notEd)').*(Tiso.*(1-phv) + dis.Ts_v1(notEd)'.*phv);
nissym        = g(:,1+2*lx+4*ln+notEd);
hospts        = g(:,1+2*lx+5*ln+notEd);
deaths        = g(:,1+2*lx+6*ln+notEd);
abs_ill       = isoasyu + isoasyv + isosymu + isosymv + nissym + hospts + deaths;%number of workers absent
abs_ill_pc    = max(0,abs_ill./data.NNs(notEd)');%proportion of workers absent
abs_ill_pcop  = abs_ill_pc.*g(:,1+notEd);%proportion of workers absent in open sectors
abs_ilop_int  = trapz(t,abs_ill_pcop);%proportion of worker-days of absence
gdpl_lbs      = abs_ilop_int.*data.obj(notEd)';
cost(7,notEd) = gdpl_lbs;

%labour demand: GDP loss due to business closures
x             = g(:,1+notEd);%x = w.^data.alp;%proportion of sector open
x_int         = trapz(t,1-x);%diff(t)'*(1-x(1:end-1,:));%proportion of sector closed (summed)
gdpl_lbd      = x_int.*data.obj(notEd)';
cost(8,notEd) = gdpl_lbd;

%% SUMMARY

c1        = [cost(3,lx+1) cost(3,lx+2) sum(cost(3,[1:lx,lx+3])) cost(3,ln)];
c2        = sum(cost(7:8,1:lx),1);
map45to10 = [1,1,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,5,6,7,7,7,7,7,6,8,8,8,9,9,9,9,10,10,10,10,10,10];
c2        = accumarray(map45to10',c2')';
c3        = sum(cost(6,:));
c         = [c1,c2,c3];

end
function [cost,ccost_t] = p2Cost(data,dis,p2,g)

t  = g(:,1);
lx = length(data.obj);
ln = lx+4;

%% VLYL

deaths    = g(end,(1+2*lx+6*ln+1):(1+2*lx+7*ln));
cost(1,:) = deaths;

lyl       = deaths.*data.lgh;
cost(2,:) = lyl;

vlyl      = lyl*data.vly;
cost(3,:) = vlyl;

deaths          = g(:,(1+2*lx+6*ln+1):(1+2*lx+7*ln));
ccost_t(:,1:ln) = deaths.*data.lgh.*data.vly;

%% VSYL

Stu              = lx+2;
students         = data.NNs(Stu);
cost(4,lx+[1,2]) = students;%student numbers not used, but storing supply costs in lx+1 and demand costs in lx+2

%student supply: learning losses due to illness, when schools are open
%assume that the isolation period is the latent period plus maximum infectious period
%assume that the symptomatic period is the same duration as the infectious period (though start/end points may be different)
%not scaled by ability for sfh because the valuation used here accounts for this
ph           = g(:,1+2*lx+8*ln+1+2*ln+Stu);
phv          = (1-dis.hv1)*dis.ph(Stu);
Ts           = g(:,1+2*lx+8*ln+1+Stu);
Tiso         = dis.Tlat + max([dis.Tay,dis.Tsr,dis.Tsh]);
isoasyu      = (g(:,1+2*lx+0*ln+Stu)./dis.Tay).*Tiso;%.*(1-(1/3));
isoasyv      = (g(:,1+2*lx+1*ln+Stu)./dis.Tay).*Tiso;%.*(1-(1/3));
isosymu      = (g(:,1+2*lx+2*ln+Stu)./Ts).*(Tiso.*(1-ph) + 1.*ph);              
isosymv      = (g(:,1+2*lx+3*ln+Stu)./dis.Ts_v1(Stu)).*(Tiso.*(1-phv) + 1.*phv);
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

% %student supply: learning loss due to illness, regardless whether schools are open or not
% %assume that ths isoasy and isosym isolation period is the latent period plus maximum infectious period: note that isoasy and isosym depend on dis.Tsr and Ts repsectively
% %assume that the nissym symptomatic period is the same duration as the infectious period (though start/end points may be different)
% %not scaled by ability for sfh because the valuation used here accounts for this
% Ts           = g(:,1+2*lx+8*ln+1+Stu);
% isoasyu      = g(:,1+2*lx+0*ln+Stu).*((dis.Tlat + max(dis.Tay,dis.Tsr))./dis.Tay);%.*(1-(1/3));
% isoasyv      = g(:,1+2*lx+1*ln+Stu).*((dis.Tlat + max(dis.Tay,dis.Tsr))./dis.Tay);%.*(1-(1/3));
% isosymu      = g(:,1+2*lx+2*ln+Stu).*((dis.Tlat + max(dis.Tay,Ts))./Ts);
% isosymv      = g(:,1+2*lx+3*ln+Stu).*((dis.Tlat + max(dis.Tay,Ts))./Ts);
% nissym       = g(:,1+2*lx+4*ln+Stu);
% hospts       = g(:,1+2*lx+5*ln+Stu);
% deaths       = g(:,1+2*lx+6*ln+Stu);
% abs          = isoasyu + isoasyv + isosymu + isosymv + nissym + hospts + deaths;%numbers of students
% absint       = trapz(t,abs)/365;
% cost(5,lx+1) = absint;
% vsyl_sts     = absint*data.vsy;
% cost(6,lx+1) = vsyl_sts;
% 
% %student demand: learning loss due to school closures, for healthy students only
% pres         = students-abs;%numbers of students
% presl        = pres.*(1-g(:,1+data.EdInd));%.*(1-(1/3));
% preslint     = trapz(t,presl)/365;%(diff(t)'*presl)/365;
% cost(5,lx+2) = preslint;
% vsyl_std     = preslint*data.vsy;
% cost(6,lx+2) = vsyl_std;

%% SGDPL

notEd = [1:(data.EdInd-1),(data.EdInd+1):lx];

%labour supply: GDP loss due to illness, when businesses are open
%assume that the isolation period is the latent period plus maximum infectious period
%assume that the symptomatic period is the same duration as the infectious period (though start/end points may be different)
%asymptomatics have ability to wfh
ph            = g(:,1+2*lx+8*ln+1+2*ln+notEd);
phv           = (1-dis.hv1).*dis.ph(notEd)';
Ts            = g(:,1+2*lx+8*ln+1+notEd);
Tiso          = dis.Tlat + max([dis.Tay,dis.Tsr,dis.Tsh]);
hw            = g(:,1+1*lx+notEd);
isoasyu       = (g(:,1+2*lx+0*ln+notEd)./dis.Tay).*Tiso.*(1-hw);
isoasyv       = (g(:,1+2*lx+1*ln+notEd)./dis.Tay).*Tiso.*(1-hw);
isosymu       = (g(:,1+2*lx+2*ln+notEd)./Ts).*(Tiso.*(1-ph) + 1.*ph);              
isosymv       = (g(:,1+2*lx+3*ln+notEd)./dis.Ts_v1(notEd)').*(Tiso.*(1-phv) + 1.*phv);
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

%consumer demand
% betamod       = g(:,1+2*lx+6*ln+1);
% conloss       = min((1-betamod).*data.hconsl(notEd)',1).*data.hcon(notEd)';
% prdloss       = (absx+1-x).*data.obj(notEd)';
% difloss       = max(0,conloss-prdloss);%to avoid double counting
% gdpl_crd      = trapz(t,difloss);
cost(9,notEd) = 0;%gdpl_crd;
%zero during ld: betamod(sum(x,2)<44)=1;

%medium-term
cost(10,notEd) = 0;

ccost_t(:,2*ln+notEd) = cumtrapz(t,abs_ill_pcop,1).*data.obj(notEd)';
ccost_t(:,3*ln+notEd) = cumtrapz(t,(1-x),1).*data.obj(notEd)';

%% IMPC

% tstart       = min(p2.Tres,data.tvec(2));%response time or late lockdown time
% w            = g(:,1+notEd);
% if sum(w(end,:))==44;    
%     tend     = data.tvec(end-1);%lifting time or simulation end time
% else
%     tend     = data.tvec(end);
%     error('Implementation Cost Error!');
% end
% betamod      = g(:,1+2*lx+6*ln+1);
% hospts       = min(sum(g(:,(1+2*lx+3*ln+1):(1+2*lx+4*ln)),2),p2.Hmax);
% vaxxed       = sum(g(:,(1+2*lx+5*ln+1):(1+2*lx+6*ln)),2);%number of people vaccinated
% units        = [max(0,tend-tstart),...
%                 data.trate*sum(data.Npop/10^5)*max(0,tend-p2.t_tit),...
%                 sum(data.Npop)*trapz(t,1-betamod),...
%                 trapz(t,hospts),...
%                 2*vaxxed(end)];
% impcost      = data.pppf*(data.impcost(:,1)' + units.*data.impcost(:,2)');
% cost(11,1:5) = impcost;
% 
% cunits                = [cumtrapz(t,(tstart<t)&(t<tend)),...
%                          data.trate*sum(data.Npop/10^5)*cumtrapz(t,(p2.t_tit<t)&(t<tend)),...
%                          sum(data.Npop)*cumtrapz(t,1-betamod),...
%                          cumtrapz(t,hospts),...
%                          2*vaxxed];
% ccost_t(:,4*ln+[1:5]) = data.pppf*(data.impcost(:,1)' + cunits.*data.impcost(:,2)');

end
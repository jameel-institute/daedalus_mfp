function output = dd_run_postprocessing(data,dis,i,p2,output,tout,yout,inext)

ln    = length(data.NNs);
lx    = length(data.obj);
adInd = 3;

%% INITIALISE OUTPUT

if isempty(fieldnames(output));
    output.Tout       = data.tvec(1);
    output.Tsout      = dis.Ts';
    output.Thout      = dis.Th';
    output.phout      = dis.ph';   
    output.Iout       = zeros(1,ln);
    output.Isaout     = zeros(1,ln);
    output.Isavout    = zeros(1,ln);
    output.Issout     = zeros(1,ln);
    output.Issvout    = zeros(1,ln);
    output.Insout     = zeros(1,ln);
    output.Hout       = zeros(1,ln);
    output.Dout       = zeros(1,ln);
    output.dDout      = zeros(1,ln);
    output.Xout       = [];
    output.hwout      = [];
    output.aaout      = 0;
    output.asout      = 0;
    output.betamodout = 1;
    output.Vout       = zeros(1,ln);
end

%% CALCULATE OUTPUT AT EACH ITERATION

D     = data.Dvec(:,:,i);
S     = yout(:,0*ln+1:1*ln);
E     = yout(:,1*ln+1:2*ln);
Ina   = yout(:,2*ln+1:3*ln);
Isa   = yout(:,3*ln+1:4*ln);
Ins   = yout(:,4*ln+1:5*ln);
Iss   = yout(:,5*ln+1:6*ln);
H     = yout(:,6*ln+1:7*ln);
Shv1  = yout(:,8*ln+1:9*ln);
Sv1   = yout(:,9*ln+1:10*ln);
Ev1   = yout(:,10*ln+1:11*ln);
Inav1 = yout(:,11*ln+1:12*ln);
Isav1 = yout(:,12*ln+1:13*ln);
Insv1 = yout(:,13*ln+1:14*ln);
Issv1 = yout(:,14*ln+1:15*ln);
Hv1   = yout(:,15*ln+1:16*ln);
Sn    = yout(:,19*ln+1:20*ln);

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

occ     = sum(Hclass,2);
Hmax    = p2.Hmax;
th      = p2.th;
th0     = max(1,1+th*((occ-Hmax)/(2*Hmax-Hmax)));
pd      = min(th0.*dis.pd',1);
Th      = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
mu      = pd./Th;
dDclass = mu.*H + mu.*Hv1;     

ddk    = max(0,10^5*sum(mu.*Hclass,2)/sum(data.Npop));
sd_fun = @(a,b,c,t,d) 1./(1 + exp(a + b.*log10(d) - c.*t));%here, t is time since response
if i==1;%strcmp(data.inp3,'No Closures')||
    betamod = ones(size(occ));
elseif any(i==data.imand);
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,p2.tmand,max(p2.dmand,ddk));
else
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,tout-p2.Tres,ddk);
end

amp       = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
ph        = amp.*dis.ph';
Ts        = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
h         = ph./Ts;
cit_sw    = double((tout >= p2.t_tit) & (i~=5));
prev_sw   = double(sum(E + Ev1 + Ina + Isa + Ins + Iss + Inav1 + Isav1 + Insv1 + Issv1 + H + Hv1, 2) < 10^-7*sum(data.Npop));
Tss_eff   = max(0,p2.iisym-dis.Tlat);
Tsa_eff   = max(0,p2.iitra-dis.Tlat);
tm_a      = 1*(1-cit_sw.*prev_sw) + (min(Tsa_eff,dis.Tay)./dis.Tay)*(cit_sw.*prev_sw);
tm_s      = 1*(1-cit_sw) + (min(Tss_eff,Ts)./Ts).*(cit_sw.*(1-prev_sw)) + (min(Tsa_eff,Ts)./Ts).*(cit_sw.*prev_sw);
I         = (dis.red*Ina+Ins) + (1-dis.trv1)*(dis.red*Inav1+Insv1) + tm_a.*dis.red.*(Isa+(1-dis.trv1)*Isav1) + tm_s.*(Iss+(1-dis.trv1)*Issv1);
NN        = data.NNs;
NN(NN==0) = 1;
foi       = dis.beta*betamod.*(D*(I'./NN))';
seedvec   = 1e-8*sum(data.Npop)*dis.Ev*data.xconf(i,data.IntlInd);%one billionth of the population
seed      = dis.beta*betamod.*(D*(seedvec./NN))';
inc       = sum(S.*(foi+seed) + Shv1.*(foi+seed) + Sv1.*(1-dis.scv1).*(foi+seed), 2);
sinc      = dis.ps*inc;
hadm      = sum(h.*Ins + h.*Iss + dis.h_v1'.*Insv1 + dis.h_v1'.*Issv1, 2);
asc_a     = 0*(1-cit_sw.*prev_sw) + min(p2.masc, max(0, p2.trate - hadm)./inc).*(cit_sw.*prev_sw);
asc_s     = 0*(1-cit_sw) + min(p2.masc, max(0, p2.trate - hadm)./sinc).*(cit_sw.*(1-prev_sw)) + min(p2.masc, max(0, p2.trate - hadm)./inc).*(cit_sw.*prev_sw);

X  = data.xconf(i,:).*ones(length(tout),lx);
hw = data.hw(i,:).*ones(length(tout),lx);

output.dmand = ddk(end);

%% STORE OUTPUT AT EACH ITERATION 

output.Tout       = [output.Tout;tout(2:end)];  
output.Tsout      = [output.Tsout;Ts(2:end,:)];  
output.Thout      = [output.Thout;Th(2:end,:)];
output.phout      = [output.phout;ph(2:end,:)];
output.Iout       = [output.Iout;Iclass(2:end,:)];
output.Isaout     = [output.Isaout;Isaclass(2:end,:)];
output.Isavout    = [output.Isavout;Isavlass(2:end,:)];
output.Issout     = [output.Issout;Issclass(2:end,:)];
output.Issvout    = [output.Issvout;Issvlass(2:end,:)];
output.Insout     = [output.Insout;Insclass(2:end,:)];
output.Hout       = [output.Hout;Hclass(2:end,:)];
output.Dout       = [output.Dout;Dclass(2:end,:)];
output.dDout      = [output.dDout;dDclass(2:end,:)]; 
output.Xout       = [output.Xout;X(1:end-1,:)];      
output.hwout      = [output.hwout;hw(1:end-1,:)];
output.aaout      = [output.aaout;asc_a(2:end)];
output.asout      = [output.asout;asc_s(2:end)];
output.betamodout = [output.betamodout;betamod(2:end)];
output.Vout       = [output.Vout;Vclass(2:end,:)];

%% FINAL OUTPUT

if inext==6;
    output.Xout  = [output.Xout;output.Xout(end,:)];
    output.hwout = [output.hwout;output.hwout(end,:)];
    output.f     = [output.Tout,output.Xout,output.hwout,...
                    output.Isaout,output.Isavout,output.Issout,output.Issvout,output.Insout,output.Hout,output.Dout,output.Vout,...
                    output.betamodout,output.Tsout,output.Thout,output.phout];
    output.g     = [output.Tout,sum(output.Iout,2),sum(output.Hout,2),sum(output.Dout,2),...
                    output.aaout,output.asout,output.betamodout,...  
                    sum(output.Vout(:,lx+1),2),sum(output.Vout(:,lx+2),2),sum(output.Vout(:,[1:lx,lx+adInd]),2),sum(output.Vout(:,lx+4),2),...
                    sum(output.Dout(:,lx+1),2),sum(output.Dout(:,lx+2),2),sum(output.Dout(:,[1:lx,lx+adInd]),2),sum(output.Dout(:,lx+4),2),...
                    sum(output.dDout(:,lx+1),2),sum(output.dDout(:,lx+2),2),sum(output.dDout(:,[1:lx,lx+adInd]),2),sum(output.dDout(:,lx+4),2)];
end

end
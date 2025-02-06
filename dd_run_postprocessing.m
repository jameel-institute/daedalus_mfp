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
    
occ   = sum(Hclass,2);
Hmax  = p2.Hmax;
th    = p2.th;
th0   = max(1,1+th*((occ-Hmax)/(2*Hmax-Hmax)));

pd  = min(th0.*dis.pd',1);
Th  = ((1-pd).*dis.Threc)+(pd.*dis.Thd);
mu  = pd./Th;
ddk = max(0,10^5*sum(mu.*Hclass,2)/sum(data.Npop));

sd_fun = @(a,b,c,t,d) 1./(1 + exp(a + b.*log10(d) - c.*t));%here, t is time since response

if i==1;%strcmp(data.inp3,'No Closures')||
    betamod = ones(size(occ));
elseif any(i==data.imand);
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,p2.tmand,max(p2.dmand,ddk));
else
    betamod = sd_fun(p2.sda,p2.sdb,p2.sdc,tout-p2.Tres,ddk);
end

S     = yout(:,0*ln+1:1*ln);
Sn    = yout(:,19*ln+1:20*ln);
Ins   = yout(:,4*ln+1:5*ln);
Iss   = yout(:,5*ln+1:6*ln);
Insv1 = yout(:,13*ln+1:14*ln);
Issv1 = yout(:,14*ln+1:15*ln);
H     = yout(:,6*ln+1:7*ln);
Hv1   = yout(:,15*ln+1:16*ln);

amp   = min((Sn+(1-dis.heff).*(S-Sn))./S,1);
ph    = amp.*dis.ph';
Ts    = ((1-ph).*dis.Tsr) + (ph.*dis.Tsh);
g3    = (1-pd)./Th;
h     = ph./Ts;
h_v1  = dis.h_v1;

Hdot   = h.*Ins      +h.*Iss        -(g3+mu).*H;
Hv1dot = h_v1'.*Insv1 +h_v1'.*Issv1   -(g3+mu).*Hv1;
occdot = sum(Hdot+Hv1dot,2);
r      = occdot./occ;

dDclass = mu.*H + mu.*Hv1;     

trate  = p2.trate;
asca   = p2.asca;
ascb   = p2.ascb;
ascc   = p2.ascc;
pcta   = p2.pcta;
pctb   = p2.pctb;
opsa   = p2.opsa;
opsb   = p2.opsb;
opc    = p2.opc;
ps     = dis.ps;
siga   = dis.siga;
sigs   = dis.sigs;
E      = yout(:,1*ln+1:2*ln);
Ev1    = yout(:,10*ln+1:11*ln);
incid  = max(0,10^5*((siga+sigs)*sum(E+Ev1,2))/sum(data.Npop));

if i~=5;
    asc_s  = 1./(1+exp(asca+ascb*log10(incid)+ascc*log10(trate)));
    propCT = 1./(1+exp(pcta+pctb*log10(incid)));
    asc_a  = propCT.*asc_s + (1-propCT)*0;
    
    asc_a  = min(asc_a,(trate./incid).*(0*(1-propCT) + (1-ps)*propCT));
    asc_s  = min(asc_s,(trate./incid).*(1*(1-propCT) + ps*propCT));
    asc_a  = max(trate/10^5*(0*(1-propCT) + (1-ps)*propCT),asc_a);
    asc_s  = max(trate/10^5*(1*(1-propCT) + ps*propCT),asc_s);    

    asc_a  = asc_a.*(tout>=p2.t_tit);    
    asc_s  = asc_s.*(tout>=p2.t_tit);    

else
    asc_a  = zeros(size(tout));
    asc_s  = zeros(size(tout));
end

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
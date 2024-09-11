clist      = dir(fullfile('*.mat'));
clist      = {clist.name};
clist(end) = [];
X          = zeros(length(clist),9);
for i = 1:length(clist);
    load(clist{i},'data');  
    X(i,1) = data.Tres;
    X(i,2) = data.t_tit;
    X(i,3) = data.trate;
    X(i,4) = data.sdl;
    X(i,5) = data.sdb;
    X(i,6) = data.Hmax;
    X(i,7) = data.t_vax;
    X(i,8) = data.arate;
    X(i,9) = data.puptake;   
end

var = X(:,9);
lim = [0,1];
figure;
histogram((var),length(clist));

p   = fitdist((var),'normal');
figure;
hold on;
sup = linspace(min(var),max(var),1000);
plot(sup,pdf(p,sup),'r');

pt = truncate(p,lim(1),lim(2));
plot(sup,pdf(pt,sup),'g');

p2.Tres  =      random(truncate(makedist('normal','mu',66.1964,'sigma',12.3309),-75,Inf));
p2.t_tit =      random(truncate(makedist('lognormal','mu',4.12346,'sigma',0.565271),-75,Inf));
p2.trate =      random(truncate(makedist('lognormal','mu',5.16249,'sigma',1.17374),0,Inf));   
p2.sdl   =      random(truncate(makedist('uniform','lower',0.1,'upper',0.4),0.1,0.4)); 
p2.sdb   =  exp(random(truncate(makedist('normal','mu',-12.8411,'sigma',10.8119),-Inf,Inf)));%%%
p2.Hmax  =      random(truncate(makedist('lognormal','mu',4.12546,'sigma',0.836187),0,Inf))*sum(data.Npop)/(10^5);
t_vax    =      random(truncate(makedist('normal','mu',450.916,'sigma',57.0053),-75,Inf));
arate    =      random(truncate(makedist('normal','mu',329.448,'sigma',112.91),0,Inf))*sum(data.Npop)/(10^5);
puptake  =      random(truncate(makedist('normal','mu',0.70855,'sigma',0.142042),0,1));
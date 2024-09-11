X   = readtable('country_data.csv','PreserveVariableNames', true);
lab = X.igroup;
klic       = find(strcmp(lab,'LIC'));
klmc       = find(strcmp(lab,'LMIC'));
kumc       = find(strcmp(lab,'UMIC'));
khic       = find(strcmp(lab,'HIC'));
var = table2array(X(:,4:24));var=var./sum(var,2);%var = X.trate;

figure; hold on;
plot((var(klic,:)'),'Color',[0 1 1]);
plot((var(klmc,:)'),'Color',[0 .5 1]);
plot((var(kumc,:)'),'Color',[0 0 .5]);
plot((var(khic,:)'),'Color',[0.41 0.16 0.38]);
%plot(cumsum(var'),'Color',[0 1 1]);

otest= [zeros(197,1),var(:,:)]';
otest= otest(:);
test = cumsum([zeros(197,1),var(:,:)]');
test = test(:);
testx = repmat([0:21]',197,1);
fun        = @(b0,b1,x) 1./(1+exp(b0+b1*x));
[fun_fit] = fit(testx,test,fun,'StartPoint',[-0.1,-0.5],'Lower',[-10,-10],'Upper',[10,10],...
                      'Robust','LAR','MaxIter',10^3,'MaxFunEvals',10^3);
dist = fun(fun_fit.b0,fun_fit.b1,0:21);
                  
figure;hold on;
plot(testx,test,'r*');
plot(0:21,dist,'b-');

figure;hold on;
plot(testx,otest,'r*');
plot(0:21,[dist(1),diff(dist)],'b-');




lim = [0,1];
figure;
histogram((var),length(clist));

p   = fitdist((var),'exponential');
figure;
hold on;
sup = linspace(min(var),max(var),1000);
plot(sup,pdf(p,sup),'r');

pt = truncate(p,lim(1),lim(2));
plot(sup,pdf(pt,sup),'g');
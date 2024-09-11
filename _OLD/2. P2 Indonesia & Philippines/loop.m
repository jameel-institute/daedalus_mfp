close all;
clear all;

countries = {'Indonesia','Philippines'};
diseases  = {'Swine Flu','Spanish Flu','Covid','SARS'};
p2levels  = {'BAD','MEDIUM','GOOD'};
d         = zeros(length(countries),length(diseases),length(p2levels));
g         = zeros(length(countries),length(diseases),length(p2levels));
dn        = zeros(length(countries),length(diseases),length(p2levels));
gn        = zeros(length(countries),length(diseases),length(p2levels));

for i=1:length(countries)
    for j=1:length(diseases)
        for k=1:length(p2levels)
            
            country   = countries{i};
            disease   = diseases{j};
            p2level   = p2levels{k};
            h         = p2Sim_RC(country,disease,p2level);
            d(i,j,k)  = h(1);
            g(i,j,k)  = h(2);
            dn(i,j,k) = h(3);
            gn(i,j,k) = h(4);
            
        end
    end    
end

cost_india = 10333;
cost_indon = cost_india*(270624/1366417);
cost_phili = cost_india*(108115/1366417);

gdp_indon = 1073115.98;
gdp_phili = 308333.80;

%% Table 5

table5 = [dn(1,1,2)   gn(1,1,2);...
          dn(2,1,2)   gn(2,1,2);...
          dn(1,2,2)   gn(1,2,2);...
          dn(2,2,2)   gn(2,2,2);...
          dn(1,3,2)   gn(1,3,2);...
          dn(2,3,2)   gn(2,3,2);...
          dn(1,4,2)   gn(1,4,2);...
          dn(2,4,2)   gn(2,4,2)];

writematrix(table5,'table5.csv');

%% Table 6

table6 = [dn(1,1,1)-dn(1,1,3),  dn(1,2,1)-dn(1,2,3),    dn(1,3,1)-dn(1,3,3),    dn(1,4,1)-dn(1,4,3);
          dn(2,1,1)-dn(2,1,3),  dn(2,2,1)-dn(2,2,3),    dn(2,3,1)-dn(2,3,3),    dn(2,4,1)-dn(2,4,3)];

writematrix(table6,'table6.csv');

%% Table 7

table7 = [g(1,1,1)-g(1,1,3),    cost_indon, g(1,1,1)-g(1,1,3)-cost_indon,   100*(g(1,1,1)-g(1,1,3)-cost_indon)/gdp_indon;...
          g(2,1,1)-g(2,1,3),    cost_phili, g(2,1,1)-g(2,1,3)-cost_phili,   100*(g(2,1,1)-g(2,1,3)-cost_phili)/gdp_phili;...
          ...
          g(1,2,1)-g(1,2,3),    cost_indon, g(1,2,1)-g(1,2,3)-cost_indon,   100*(g(1,2,1)-g(1,2,3)-cost_indon)/gdp_indon;...
          g(2,2,1)-g(2,2,3),    cost_phili, g(2,2,1)-g(2,2,3)-cost_phili,   100*(g(2,2,1)-g(2,2,3)-cost_phili)/gdp_phili;...
          ...
          g(1,3,1)-g(1,3,3),    cost_indon, g(1,3,1)-g(1,3,3)-cost_indon,   100*(g(1,3,1)-g(1,3,3)-cost_indon)/gdp_indon;...
          g(2,3,1)-g(2,3,3),    cost_phili, g(2,3,1)-g(2,3,3)-cost_phili,   100*(g(2,3,1)-g(2,3,3)-cost_phili)/gdp_phili;...
          ...
          g(1,4,1)-g(1,4,3),    cost_indon, g(1,4,1)-g(1,4,3)-cost_indon,   100*(g(1,4,1)-g(1,4,3)-cost_indon)/gdp_indon;...
          g(2,4,1)-g(2,4,3),    cost_phili, g(2,4,1)-g(2,4,3)-cost_phili,   100*(g(2,4,1)-g(2,4,3)-cost_phili)/gdp_phili];

writematrix(table7,'table7.csv');

%% Table 9 

table9a = [(g(1,1,1)-g(1,1,3))/cost_indon,(g(1,2,1)-g(1,2,3))/cost_indon,(g(1,3,1)-g(1,3,3))/cost_indon,(g(1,4,1)-g(1,4,3))/cost_indon;...   
           (g(2,1,1)-g(2,1,3))/cost_phili,(g(2,2,1)-g(2,2,3))/cost_phili,(g(2,3,1)-g(2,3,3))/cost_phili,(g(2,4,1)-g(2,4,3))/cost_phili];  

writematrix(table9a,'table9a.csv');

table9b = 10^6*1./table9a;

writematrix(table9b,'table9b.csv');

%%

f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.2 0.2 0.70 0.70];
hold on;

inds=[1:2];
data=...
[6.95 210.05 289.32 396.33; ...
 4.49 129.13 164.40 261.54];
scatter(inds,data(:,1),50,'blue','filled','^');
scatter(inds,data(:,2),50,'red','filled','o');
scatter(inds,data(:,3),50,'magenta','filled','s');
scatter(inds,data(:,4),50,'black','filled','d');
set(gca,'YScale','log');

axis square;
box on;
grid on;
grid minor;
xticks(inds);
xticklabels({'Indonesia','Philippines'});
yticks([1,10,100,1000]);
xlim([0.5 2.5]);
ylim([1 5000]);
xlabel('Country');
ylabel('GDP Loss Averted/\$ P2 Cost');
vec_pos=get(get(gca,'xlabel'),'Position');
set(get(gca,'xlabel'),'Position',vec_pos+[0 -0.04 0]);        
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.15 0 0]);
legend('Swine flu','Spanish flu','COVID','SARS','position',[0.696 0.76 0.1 0.1]);
vec_pos=get(get(gca,'title'),'Position');
set(get(gca,'title'),'Position',vec_pos+[-3 5000 0]);
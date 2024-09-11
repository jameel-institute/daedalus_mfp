f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.2 0.2 0.70 0.70];
hold on;

inds=[1:4];
data=...
[-9.555240297	81.43061901	621.1058732	222.6151507;
 -9.981029433	64.89999668	646.9976719	200.9268413;
 -8.088061141	85.82584601	467.3571072	216.415833;
 -2.859946054	6.82534217	147.0433352	99.06805129];
scatter(inds,data(:,2),50,'red','filled','o');
scatter(inds,data(:,3),50,'magenta','filled','s');
scatter(inds,data(:,4),50,'black','filled','d');

axis square;
box on;
grid on;
grid minor;
xticks(inds);
xticklabels({'USA','UK','China','India'});
xlim([0.5 4.5]);
ylim([0 700]);
xlabel('Country');
ylabel('Deaths Averted/\$1m P2 Cost');
vec_pos=get(get(gca,'xlabel'),'Position');
set(get(gca,'xlabel'),'Position',vec_pos+[0 -35 0]);
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.15 0 0]);
legend('Spanish flu','COVID','SARS','position',[0.696 0.78 0.1 0.1]);
title('(a)');
vec_pos=get(get(gca,'title'),'Position');
set(get(gca,'title'),'Position',vec_pos+[-3 50 0]);

%%

f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.2 0.2 0.70 0.70];
hold on;

inds=[1:4];
data1=...
[-28.27803995	414.6041088	1702.835294	967.6062552;
 -18.80064173	249.084372	1095.62268	610.7795473;
 -5.825620811	114.5712097	295.2935345	244.6520591;
 -0.89551349	4.318438873	33.59916895	41.68418664];
data2=...
[-45.24486392	663.3665741	2724.53647	1548.170008;
 -30.08102676	398.5349952	1752.996287	977.2472756;
 -9.320993298	183.3139355	472.4696551	391.4432946;
 -1.432821584	6.909502196	53.75867031	66.69469862];
scatter(inds+0.15,data1(:,2),50,'red','filled','o');
scatter(inds-0.15,data1(:,3),50,'magenta','filled','s');
scatter(inds,data1(:,4),50,'black','filled','d');
scatter(inds+0.15,data2(:,2),50,'red','filled','o');
scatter(inds-0.15,data2(:,3),50,'magenta','filled','s');
scatter(inds,data2(:,4),50,'black','filled','d');
for i=1:4
    plot((inds(i)+0.15)*[1,1],[data1(i,2),data2(i,2)],'-r','linewidth',2);
    plot((inds(i)-0.15)*[1,1],[data1(i,3),data2(i,3)],'-m','linewidth',2);
    plot((inds(i))*[1,1],[data1(i,4),data2(i,4)],'-k','linewidth',2);
end
set(gca,'YScale','log');

axis square;
box on;
grid on;
grid minor;
xticks(inds);
xticklabels({'USA','UK','China','India'});
yticks([1,10,100,1000,5000]);
xlim([0.5 4.5]);
ylim([0.25 5000]);
xlabel('Country');
ylabel('VYLL Averted/\$ P2 Cost');
vec_pos=get(get(gca,'xlabel'),'Position');
set(get(gca,'xlabel'),'Position',vec_pos+[0 -0.04 0]);                   
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.15 0 0]);
legend('Spanish flu','COVID','SARS','position',[0.696 0.78 0.1 0.1]);
title('(b)');
vec_pos=get(get(gca,'title'),'Position');
set(get(gca,'title'),'Position',vec_pos+[-3 5000 0]);                   

%%

f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.2 0.2 0.70 0.70];
hold on;

inds=[1:4];
data=...
[23.29158001	42.98277116	1102.426348	82.78393378;
 15.66602677	29.09399087	788.9300698	55.49533361;
 3.729752905	6.820492737	164.0833051	13.26676625;
 0.28741853	    0.760128659	9.68221412	1.411120429];
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
xticklabels({'USA','UK','China','India'});
yticks([1,10,100,1000,5000]);
xlim([0.5 4.5]);
ylim([0.25 5000]);
xlabel('Country');
ylabel('GDP Loss Averted/\$ P2 Cost');
vec_pos=get(get(gca,'xlabel'),'Position');
set(get(gca,'xlabel'),'Position',vec_pos+[0 -0.04 0]);        
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.15 0 0]);
legend('Swine flu','Spanish flu','COVID','SARS','position',[0.696 0.76 0.1 0.1]);
title('(c)');
vec_pos=get(get(gca,'title'),'Position');
set(get(gca,'title'),'Position',vec_pos+[-3 5000 0]);
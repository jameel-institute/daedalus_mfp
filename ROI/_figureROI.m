%close all;
%clear all;

countries  = {'Ethiopia',...
             'India','Philippines','Indonesia',...
             'South Africa','Brazil','China',...
             'United Kingdom','Australia','United States'};
diseases   = {'Flu2009','Flu1957','Flu1918',...
             'CovidWT','CovidOM','CovidDE',...
             'SARS'};
strategies = {'Unmitigated','School Closures','Economic Closures','Elimination'};

%%

Npop = zeros(length(countries),1);
Avag = zeros(length(countries),1);
GDP  = zeros(length(countries),1);
EPVc = zeros(length(countries),5);

for i = 1:length(countries);
    
    load(strcat(countries{i},'.mat'),'data');
    Npop(i)   = sum(data.Npop);
    Avag(i)   = dot(data.Npop,[2:5:102]')/Npop(i);
    GDP(i)    = sum(365*data.obj);
    costs     = data.pppf*data.prepcost;                 
    startu    = costs(:,1:2:5);
    annual    = costs(:,2:2:6);
    startu    = sum(startu(:,3:4-1),2);%levels
    annual    = sum(annual(:,3:4-1),2);%levels
    n         = 10;%check same below
    r         = 0.03;             
    EPVc(i,:) = (startu + annual*sum(1./(1+r).^[0:n-1]))';
    
end

%%

EPVg = zeros(length(countries),length(diseases),length(strategies));
ROIl = zeros(length(countries),length(diseases),length(strategies));
ROI  = zeros(length(countries),length(diseases),length(strategies));
ROIu = zeros(length(countries),length(diseases),length(strategies));

for i = 1;%[1,6,10];%:length(countries);
    for j = 6;%:length(diseases);
        for k = 1:length(strategies);

            [ROIl(i,j,k),ROI(i,j,k),ROIu(i,j,k),EPVg(i,j,k)] = p2ROI(countries{i},diseases{j},strategies{k},'LEVEL3','LEVEL4',0.0041,0.0048,0.0055,n);

            disp(i);
            
        end
    end
end

%%

f  = figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

ax = gca;
ax.Position = [0.08 0.28 0.88 0.625];
hold on;

labs = categorical(countries);
labs = reordercats(labs,[4 5 7 6 8 2 3 9 1 10]);
y    = Avag;
b    = bar(labs,y,'stacked','FaceColor','flat');

%ylim([0 ylmt]);
%axis square;
xtickangle(45);
grid on;
grid minor;
box on;
ylabel('Average Age');
set(gca,'FontSize',fs);

b.CData(1,:)    = repmat([0 1 1],1,1);
b.CData(2:4,:)  = repmat([0 .5 1],3,1);
b.CData(5:7,:)  = repmat([0 0 .5],3,1);
b.CData(8:10,:) = repmat([0.41 0.16 0.38],3,1);
%legend('Response Time','Testing \& Tracing','Social/Physical Distancing','Hospital Capacity','Vaccination','location','northwest');

%%

f  = figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

ax = gca;
ax.Position = [0.08 0.28 0.88 0.625];
hold on;

labs = categorical(countries);
labs = reordercats(labs,[4 5 7 6 8 2 3 9 1 10]);
y    = 10^6*EPVc./Npop;
b    = bar(labs,y,'stacked','FaceColor','flat');

%ylim([0 ylmt]);
%axis square;
xtickangle(45);
grid on;
grid minor;
box on;
ylabel('Level 3-4 Preparedness Costs pc (\$)');
set(gca,'FontSize',fs);

b(1).CData = repmat([0.0000    0.4470    0.7410],length(countries),1);
b(2).CData = repmat([0.8500    0.3250    0.0980],length(countries),1);
b(3).CData = repmat([0.9290    0.6940    0.1250],length(countries),1);
b(4).CData = repmat([0.4940    0.1840    0.5560],length(countries),1);
b(5).CData = repmat([0.4660    0.6740    0.1880],length(countries),1);
legend('Response Time','Testing \& Tracing','Social/Physical Distancing','Hospital Capacity','Vaccination','location','northwest');

%%

f  = figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

ax = gca;
ax.Position = [0.08 0.2 0.88 0.70];
hold on;

xmaj  = repelem(1:length(diseases),length(strategies));
xmin  = repmat([-0.3 -0.1 0.1 0.3],1,length(diseases));
x     = xmaj + xmin;
y     = 10^6*EPVg./Npop;
z     = zeros(length(countries),length(diseases)*length(strategies));
for i = 1:length(countries);
    zz     = squeeze(y(i,:,:));
    z(i,:) = reshape(zz',[],1)';
end
style = {'s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o'};
color = jet(13);
color = color([1,3,5,7:13],:);
h1    = plot(x,z(1,:),'o','MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:));
h2    = plot(x,z(2,:),'o','MarkerEdgeColor',color(2,:),'MarkerFaceColor',color(2,:));
h3    = plot(x,z(3,:),'o','MarkerEdgeColor',color(3,:),'MarkerFaceColor',color(3,:));
h4    = plot(x,z(4,:),'o','MarkerEdgeColor',color(4,:),'MarkerFaceColor',color(4,:));
h5    = plot(x,z(5,:),'o','MarkerEdgeColor',color(5,:),'MarkerFaceColor',color(5,:));
h6    = plot(x,z(6,:),'o','MarkerEdgeColor',color(6,:),'MarkerFaceColor',color(6,:));
h7    = plot(x,z(7,:),'o','MarkerEdgeColor',color(7,:),'MarkerFaceColor',color(7,:));
h8    = plot(x,z(8,:),'o','MarkerEdgeColor',color(8,:),'MarkerFaceColor',color(8,:));
h9    = plot(x,z(9,:),'o','MarkerEdgeColor',color(9,:),'MarkerFaceColor',color(9,:));
h10   = plot(x,z(10,:),'o','MarkerEdgeColor',color(10,:),'MarkerFaceColor',color(10,:));

set(gca,'YScale','log');
xlim([0 9]);
%ylim([]);
%axis square;
xticks(1:7);
grid on;
grid minor;
box on;
xticklabels(diseases);
xtickangle(45);
%xlabel('Time');
ylabel('Level 3-4 Preparedness Gains pc (\$)');
set(gca,'FontSize',fs);

legend(countries);
%title(countries(i));

%%

f  = figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

ax = gca;
ax.Position = [0.08 0.2 0.88 0.70];
hold on;

x     = [1:7];    
xp    = [-0.3 -0.1 0.1 0.3];
style = {'s','d','^','o'};
%color = parula(10);
color = [0.00 1.00 0.00;
         0.93 0.69 0.13;
         1.00 0.50 0.00;
         1.00 0.00 1.00;
         1.00 0.00 0.00; 
         0.00 0.00 1.00;
         0.00 0.00 0.00];

for i = 1;%1:length(countries);
    for j = 1:length(diseases);
        
        y  = ROI(i,j,:);
        en = ROI(i,j,:) -ROIl(i,j,:);
        ep = ROIu(i,j,:)-ROI(i,j,:);
        
        for k = 1:length(strategies);
            h = errorbar(x(j)+xp(k),y(k),en(k),ep(k),style{k},'Color',color(j,:),'MarkerFaceColor',color(j,:),'MarkerSize',7);
        end
        
    end
end
plot(0:8,7*ones(1,9),'-','LineWidth',1,'Color','red');

%set(gca,'YScale','log');
%ylim([10^1 1.4*10^5]);
%axis square;
xticks(1:7);
grid on;
grid minor;
box on;
xticklabels(diseases);
xtickangle(45);
ax.YAxis.Exponent = 0;
%xlabel('Time');
ylabel('Level 3-4 ROI (\%)');
set(gca,'FontSize',fs);

title(countries(i));

%%

f  = figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

ax = gca;
ax.Position = [0.08 0.26 0.88 0.65];
hold on;

x     = repelem(1:length(countries),length(strategies));
y     = 10^6*EPVg/1.07./Npop;
z     = zeros(length(diseases),length(countries)*length(strategies));
for i = 1:length(diseases);
    zz     = squeeze(y(:,i,:));
    z(i,:) = reshape(zz',[],1)';
end
% style = {'s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o','s','d','^','o'};
color = [0.00 1.00 0.00;
         0.93 0.69 0.13;
         1.00 0.50 0.00;
         1.00 0.00 1.00;
         1.00 0.00 0.00; 
         0.00 0.00 1.00;
         0.00 0.00 0.00];
h1    = plot(x,z(1,:),'o','MarkerEdgeColor',color(1,:),'MarkerFaceColor',color(1,:));
h2    = plot(x,z(2,:),'o','MarkerEdgeColor',color(2,:),'MarkerFaceColor',color(2,:));
h3    = plot(x,z(3,:),'o','MarkerEdgeColor',color(3,:),'MarkerFaceColor',color(3,:));
h4    = plot(x,z(4,:),'o','MarkerEdgeColor',color(4,:),'MarkerFaceColor',color(4,:));
h5    = plot(x,z(5,:),'o','MarkerEdgeColor',color(5,:),'MarkerFaceColor',color(5,:));
h6    = plot(x,z(6,:),'o','MarkerEdgeColor',color(6,:),'MarkerFaceColor',color(6,:));
h7    = plot(x,z(7,:),'o','MarkerEdgeColor',color(7,:),'MarkerFaceColor',color(7,:));
 
set(gca,'YScale','log');
xlim([0 12]);
%ylim([]);
%axis square;
xticks(1:10);
grid on;
grid minor;
box on;
xticklabels(countries);
xtickangle(45);
%xlabel('Time');
ylabel('Level 3-4 Preparedness Budget pc (\$)');
set(gca,'FontSize',fs);

legend(diseases);
%title(countries(i));
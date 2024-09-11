close all;
clear all;

addpath('../../../../../');
load('IN36.mat','data');

%%

% LASTN = maxNumCompThreads(1);

% par1=linspace(105,1210,60);
% par2=linspace(0.60,0.95,60);
% 
% deat=zeros(length(par1),length(par2));
% gdpl=zeros(length(par1),length(par2));
% 
% i=1;
% for inp1=par1   
%     j=1;
%     for inp2=par2
%     
%     [~,h2]=heSwitchSim(inp1,inp2);
%     deat(i,j)=h2(1);
%     gdpl(i,j)=h2(2);
%     
%     j=j+1;
%     end    
%     i=i+1
% end
% file=[par1;par2;deat;gdpl];
% dlmwrite('imm_rsk.txt',file);
% %gdpl=100*gdpl/sum(2*12*data.obj/1000);

%%

file=dlmread('imm_rsk.txt');
par1=file(1,:);
par2=file(2,:);
deat=file(3:2+length(par1),:);
gdpl=file(3+length(par1):end,:);
gdpl=100*gdpl/sum(2*12*data.obj/1000);

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(par1,100*par2,deat');
%plot3(682.5457,75.42,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(650.1573,71.48,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(152.3081,88.62,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(105.1017,74.53,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([par1(1) par1(end)]);
ylim([100*par2(1) 100*par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Administration Rate (per 100k/day)');
ylabel('Uptake (\%)');
colormap(customcolormap_preset('red-yellow-green'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{Deaths (per 100k)}');

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(par1,100*par2,gdpl');
%plot3(682.5457,75.42,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(650.1573,71.48,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(152.3081,88.62,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(105.1017,74.53,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([par1(1) par1(end)]);
ylim([100*par2(1) 100*par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Administration Rate (per 100k/day)');
ylabel('Uptake (\%)');
colormap(customcolormap_preset('red-yellow-blue'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{GDP Loss (\% per biennium)}');

%%

% LASTN = maxNumCompThreads(1);

% par1=linspace(0.60,0.40,60);
% par2=linspace(12,252,60);
% 
% deat=zeros(length(par1),length(par2));
% gdpl=zeros(length(par1),length(par2));
% 
% i=1;
% for inp1=par1   
%     j=1;
%     for inp2=par2
%     
%     [~,h2]=heSwitchSim(inp1,inp2);
%     deat(i,j)=h2(1);
%     gdpl(i,j)=h2(2);
%     
%     j=j+1;
%     end    
%     i=i+1
% end
% file=[par1;par2;deat;gdpl];
% dlmwrite('npl_mcm.txt',file);
% %gdpl=100*gdpl/sum(2*12*data.obj/1000);

%%

file=dlmread('npl_mcm.txt');
par1=file(1,:);bindex=@(del) (1-del)*1000/6;
par2=file(2,:);
deat=file(3:2+length(par1),:);
gdpl=file(3+length(par1):end,:);
gdpl=100*gdpl/sum(2*12*data.obj/1000);

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(bindex(par1),par2,deat');
%plot3(bindex(0.5801),93.6250,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(bindex(0.54),33.1408,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(bindex(0.5267),66.4339,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(bindex(0.4223),22.9357,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([bindex(par1(1)) bindex(par1(end))]);
ylim([par2(1) par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Adherence to NPIs');
ylabel('Hospital Capacity (per 100k)');
colormap(customcolormap_preset('red-yellow-green'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{Deaths (per 100k)}');

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(bindex(par1),par2,gdpl');
%plot3(bindex(0.5801),93.6250,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(bindex(0.54),33.1408,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(bindex(0.5267),66.4339,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(bindex(0.4223),22.9357,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([bindex(par1(1)) bindex(par1(end))]);
ylim([par2(1) par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Adherence to NPIs');
ylabel('Hospital Capacity (per 100k)');
colormap(customcolormap_preset('red-yellow-blue'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{GDP Loss (\% per biennium)}');

%%

% LASTN = maxNumCompThreads(1);

% par1=linspace(-10,30,60);
% par2=linspace(0,0.025,60);
% 
% deat=zeros(length(par1),length(par2));
% gdpl=zeros(length(par1),length(par2));
% 
% i=1;
% for inp1=par1   
%     j=1;
%     for inp2=par2
%     
%     [~,h2]=heSwitchSim(inp1,inp2);
%     deat(i,j)=h2(1);
%     gdpl(i,j)=h2(2);
%     
%     j=j+1;
%     end    
%     i=i+1
% end
% file=[par1;par2;deat;gdpl];
% dlmwrite('ero_sur.txt',file);
% %gdpl=100*gdpl/sum(2*12*data.obj/1000);

%%

file=dlmread('ero_sur.txt');
par1=file(1,:);
par2=file(2,:);
deat=file(3:2+length(par1),:);
gdpl=file(3+length(par1):end,:);
gdpl=100*gdpl/sum(2*12*data.obj/1000);

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(par1,100*par2,deat');
%plot3(0,0.38,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(0,0.68,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(0,0.03,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(0,0.06,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([par1(1) par1(end)]);
ylim([100*par2(1) 100*par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Government Response Speed-Up (days)');
ylabel('Cases Self-Isolating (\% per day)');
colormap(customcolormap_preset('red-yellow-green'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{Deaths (per 100k)}');

f=figure('Units','centimeters','Position',[0 0 12.5 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
hold on;
C=surf(par1,100*par2,gdpl');
%plot3(0,0.38,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(0,0.68,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
%plot3(0,0.03,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
plot3(0,0.06,100000,'o','MarkerSize',7,'LineWidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
xlim([par1(1) par1(end)]);
ylim([100*par2(1) 100*par2(end)]);
view([0,90]);
axis square;
box on;
grid on;
%xticks([]);    
%yticks([]);
xlabel('Government Response Speed-Up (days)');
ylabel('Cases Self-Isolating (\% per day)');
colormap(customcolormap_preset('red-yellow-blue'));
shading interp;
%caxis([-0.6,0.6]);
c=colorbar;
%c.XTick=[-0.6:0.2:0.6];
title('\textbf{GDP Loss (\% per biennium)}');
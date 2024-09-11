% close all;
% clear all;

format long;

inp1 = 'United Kingdom';
inp2 = 'Flu1918';
inp4 = 'LEVEL1';

load(strcat(inp1,'.mat'),'data');
Npop = data.Npop;
Npop = [Npop(1:16);sum(Npop(17:end))];
if strcmp(inp2,'Flu2009');
    dis = p2Params_Flu2009;
elseif strcmp(inp2,'Flu1957');
    dis = p2Params_Flu1957;
elseif strcmp(inp2,'Flu1918');
    dis = p2Params_Flu1918;
elseif strcmp(inp2,'CovidWT');
    dis = p2Params_CovidWT;    
elseif strcmp(inp2,'CovidOM');
    dis = p2Params_CovidOM;    
elseif strcmp(inp2,'CovidDE');
    dis = p2Params_CovidDE;    
elseif strcmp(inp2,'SARS');
    dis = p2Params_SARS;
else
    error('Unknown Disease!');
end
R0   = dis.R0;
ps   = dis.ps;
ihr  = dis.ihr;
ifr  = dis.ifr;

tmod = linspace(1,6,32*2)./R0;
smod = zeros(1,32*2);
sfun = @(s) logspace(-3,-1,length(smod)) - Npop'/sum(Npop)*min(ifr'.*s,min(dis.ihr'.*s,min(dis.ps'.*s,1)));
smod = fsolve(sfun,ones(size(smod)));
jmax = length(smod);

%%

eccla = zeros(length(tmod),length(smod));
ecclb = zeros(length(tmod),length(smod));
ecclc = zeros(length(tmod),length(smod));
eccld = zeros(length(tmod),length(smod));
elima = zeros(length(tmod),length(smod));
elimb = zeros(length(tmod),length(smod));
elimc = zeros(length(tmod),length(smod));
elimd = zeros(length(tmod),length(smod));

parfor i = 1:length(tmod);   
    for j = 1:jmax;   

        inp5 = tmod(i);
        inp6 = smod(j);

        try
            eccl = p2Sim(inp1,inp2,'Economic Closures',inp4,inp5,inp6);
            
            eccla(i,j) = eccl(1);
            ecclb(i,j) = eccl(2);
            ecclc(i,j) = eccl(3);
            eccld(i,j) = eccl(4);
        catch
            disp('Economic Closures');
            disp(inp5);
            disp(inp6);
        end
        try
            elim = p2Sim(inp1,inp2,'Elimination',inp4,inp5,inp6);
            
            elima(i,j) = elim(1);
            elimb(i,j) = elim(2);
            elimc(i,j) = elim(3);
            elimd(i,j) = elim(4);
        catch
            disp('Elimination');
            disp(inp5);
            disp(inp6);
        end
    
    end    
end

file = [tmod;smod;eccla;ecclb;ecclc;eccld;elima;elimb;elimc;elimd];
dlmwrite(strcat(inp4,'.txt'),file);

%%

% file = dlmread(strcat(inp4,'.txt'));
% %file = dlmread('LEVEL4.txt');
% ax1  = R0*file(1,:);
% ax2  = 100*Npop'/sum(Npop)*min(ifr'.*file(2,:),ps);
% eccl = file(2+0*length(ax1)+[1:length(ax1)],:);
% %elim = file(2+length(ax1)+[1:length(ax1)],:);
% 
% f  = figure('Units','centimeters','Position',[0 0 10 10]);
% set(f,'defaulttextInterpreter','latex');
% set(f,'defaultAxesTickLabelInterpreter','latex');
% set(f,'defaultLegendInterpreter','latex');
% set(f,'defaultColorbarTickLabelInterpreter','latex');
% set(f,'DefaultAxesFontSize',12);
% fs = 12;
% 
% ax = gca;
% ax.Position = [0.15 0.15 0.80 0.66];
% hold on;
% 
% h1 = surf(ax1,ax2,(eccl./1)');
% %h2 = contour3(ax1,ax2,(eccl./elim)',[1 1],'-k','LineWidth',2);
% 
% view([0,90]);
% set(gca,'YScale','log');
% axis([1 6 0.1 10]);
% %axis square;
% grid on;
% box on;
% %xticks([]);    
% yticks([0.1,0.3,1,3,10]);
% %set(gca,'xticklabels',{'Year 1 - Jan','Jul',...
% %                       'Year 2 - Jan','Jul',...
% %                       'Year 3 - Jan','Jul'});       
% set(gca,'yticklabels',{'0.1','0.3','1','3','10'});
% xlabel('$R_0$');
% ylabel('IFR (\%)');
% set(gca,'FontSize',fs);
% 
% colormap(customcolormap_preset('red-yellow-green'));
% %colormap(flipud(bone));
% shading interp;
% %caxis([10^1,5*10^6]);
% set(gca,'ColorScale','log');
% c = colorbar;
% %c.XTick=[-0.6:0.2:0.6];
% title('\textbf{Socio-Economic Cost}');
% vec_pos = get(get(gca,'title'),'Position');
% set(get(gca,'title'),'Position',vec_pos+[0 5 0]);
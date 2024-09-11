%% VIOLINS 1

f  = figure('Units','centimeters','Position',[0 0 30 35]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918',...
               'Covid Omicron','Covid Wildtype','Covid Delta',...
               'SARS'};
% diseases   = {'Influenza 2009','Influenza 1918',...
%               'Covid Omicron','Covid Wildtype',...
%               'SARS'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
% strategies = {'No Closures','Economic Closures','Elimination'};
%ylimmax    = [1,2,12,8,8,24,8];
ylimmax    = [250,250,1200,350,450,1000,600];
% ylimmax    = [250,1200,350,450,600];
sumstats   = zeros(length(diseases),length(strategies),3);

for j = 1:length(diseases);
    for k = 1:length(strategies);
    
    inp2 = diseases(j);
    inp3 = strategies(k);
    T1   = readtable(strcat('LLMIC_',string(inp2),'_',string(inp3),'.csv'));
    T2   = readtable(strcat('UMIC_',string(inp2),'_',string(inp3),'.csv'));
    T3   = readtable(strcat('HIC_',string(inp2),'_',string(inp3),'.csv'));

%     y1   = T1.SEC/10^6;
%     y2   = T2.SEC/10^6;
%     y3   = T3.SEC/10^6;
    gdp1 = 365*sum(table2array(T1(:,331:375)),2);
    gdp2 = 365*sum(table2array(T2(:,331:375)),2);
    gdp3 = 365*sum(table2array(T3(:,331:375)),2);    
    y1   = 100*T1.SEC./gdp1;
    y2   = 100*T2.SEC./gdp2;
    y3   = 100*T3.SEC./gdp3;
    
%     Q11      = quantile(y1,0.25);
%     Q12      = quantile(y2,0.25);
%     Q13      = quantile(y3,0.25);
%     Q31      = quantile(y1,0.75);
%     Q32      = quantile(y2,0.75);
%     Q33      = quantile(y3,0.75);
%     out1     = find(y1>Q31+1.5*(Q31-Q11));
%     out2     = find(y2>Q32+1.5*(Q32-Q12));
%     out3     = find(y3>Q33+1.5*(Q33-Q13));
%     y1(out1) = NaN;
%     y2(out2) = NaN;
%     y3(out3) = NaN;
    
    %yran = linspace(0,max([y1;y2;y3]),1000);
    yran1 = linspace(0,max(y1),1000);
    yran2 = linspace(0,max(y2),1000);
    yran3 = linspace(0,max(y3),1000);
    pdf1  = ksdensity(y1,yran1,'Support','positive');%,'Bandwidth',0.1);
    pdf2  = ksdensity(y2,yran2,'Support','positive');%,'Bandwidth',0.1);
    pdf3  = ksdensity(y3,yran3,'Support','positive');%,'Bandwidth',0.1);
    
    ymax = max([pdf1,pdf2,pdf3]);
    pdf1 = pdf1/ymax/2;
    pdf2 = pdf2/ymax/2;
    pdf3 = pdf3/ymax/2;
    
    pdf1(pdf1<0.001) = 0;
    pdf2(pdf2<0.001) = 0;
    pdf3(pdf3<0.001) = 0;
    pdp1             = pdf1;
    pdp2             = pdf2;
    pdp3             = pdf3;
    pdp1(pdp1==0&circshift(pdp1,1)==0&circshift(pdp1,-1)==0) = NaN;
    pdp2(pdp2==0&circshift(pdp2,1)==0&circshift(pdp2,-1)==0) = NaN;
    pdp3(pdp3==0&circshift(pdp3,1)==0&circshift(pdp3,-1)==0) = NaN;
    
    s = subplot(length(diseases),length(strategies),k+(j-1)*length(strategies));
    d = 0.2;
    i = 0.0125+d*(0.975*1/4)+(k-1)*(0.975*1/4);
    w = (1-d)*(0.975*1/4);
    u = (d/2)*(1/7)+(7-j)*(1/7);
    h = (1-d)*(1/7);
    set(s,'position',[i u w h]);
    hold on;
    
    xval = 0.5;
    yran = yran1;
    ydat = y1;
    pdf  = pdf1;
    pdp  = pdp1;
    col  = [1 1 0];
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    xval = 1.5;
    yran = yran2;
    ydat = y2;
    pdf  = pdf2;
    pdp  = pdp2;
    col  = [1 .5 0];
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    xval = 2.5;
    yran = yran3;
    ydat = y3;
    pdf  = pdf3;
    pdp  = pdp3;
    col  = [0 0 1];
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    %plot([xval-0.1 xval+0.1],[quantile(ydata,0.05),quantile(ydata,0.05)],'Linewidth',0.5,'Color','black');
    %plot([xval xval],[quantile(ydata,0.05),quantile(ydata,0.95)],'Linewidth',0.5,'LineStyle','-','Color','black');
    %plot([xval-0.1 xval+0.1],[quantile(ydata,0.95),quantile(ydata,0.95)],'Linewidth',0.5,'Color','black');
    %plot(linspace(0.5,4.5,100),40*ones(1,100),'Linewidth',0.5,'Color','black');
    
    xlim([0 3]);
    ylm = ylimmax(j);
    ylim([0 ylm]);
    xticks([0.5:1:2.5]);
    yticks([0 ylm/2 ylm]);
    grid on;
    grid minor;
    box on;
    xticklabels({''});
    if (j==1 && k==2);
        %ylabel('Societal Cost (\$, trillion)');
        ylabel('Societal Cost (\% of GDP)');
        vec_pos = get(get(gca,'ylabel'),'Position');
        set(get(gca,'ylabel'),'Position',vec_pos + [-0.125 0 0]);
    end
    set(gca,'FontSize',fs);   
    
    %sumstats(j,k,:) = 1000*[quantile(y1,0.50),quantile(y2,0.50),quantile(y3,0.50)]; 
    sumstats(j,k,:) = [mean(y1(~isnan(y1))),mean(y2(~isnan(y2))),mean(y3(~isnan(y3)))]; 
                     
    end
end

sumstatt = [];
for j = 1:length(diseases);
    for k = 1:length(strategies);
        sumstatt = [sumstatt;squeeze(sumstats(j,k,:))'];
    end
end
%writematrix(sumstatt,'sumstatt_usd.csv');
writematrix(sumstatt,'sumstatt_pc.csv');

%% SCATTERPLOTS

f  = figure('Units','centimeters','Position',[0 0 36 27]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

diseases   = {'Influenza 2009','Influenza 1957','Influenza 1918',...
              'Covid Omicron','Covid Wildtype','Covid Delta',...
              'SARS'};
strategies = {'No Closures','School Closures','Economic Closures','Elimination'};
ylimmax    = [1,2,12,8,8,24,8];

for j = 1:length(diseases);
    for k = 1:length(strategies);
    
    inp2 = diseases(j);
    inp3 = strategies(k);
    T1   = readtable(strcat('LLMIC_',string(inp2),'_',string(inp3),'.csv'));
    T2   = readtable(strcat('UMIC_',string(inp2),'_',string(inp3),'.csv'));
    T3   = readtable(strcat('HIC_',string(inp2),'_',string(inp3),'.csv'));
    
    gdp1 = 10^3*365*sum(table2array(T1(:,331:375)),2)/(50*10^6);
    gdp2 = 10^3*365*sum(table2array(T2(:,331:375)),2)/(50*10^6);
    gdp3 = 10^3*365*sum(table2array(T3(:,331:375)),2)/(50*10^6); 
    gdpc = [gdp1;gdp2;gdp3];
    y1   = T1.SEC/10^6;
    y2   = T2.SEC/10^6;
    y3   = T3.SEC/10^6;
    y    = [y1;y2;y3];
    
    s = subplot(length(strategies),length(diseases),j+(k-1)*length(diseases));
    d = 0.2;
    i = 0.0125+d*(0.975*1/7)+(j-1)*(0.975*1/7);
    w = (1-d)*(0.975*1/7);
    u = 0.03+(d/2)*(0.97*1/4)+(4-k)*(0.97*1/4);
    h = (1-d)*(0.97*1/4);
    set(s,'position',[i u w h]);
    hold on;
    
    scatter(gdp3,y3,[],[0 0 1],'.');
    scatter(gdp2,y2,[],[1 .5 0],'.');
    scatter(gdp1,y1,[],[0.9290,0.6940,0.1250],'.');
    
    xlim([0 150]);
    ylm = ylimmax(j);
    ylim([0 ylm]);
    xticks([0:50:150]);
    yticks([0 ylm/2 ylm]);
    grid on;
    grid minor;
    box on;
    %xticklabels({''});
    if (j==4 && k==4);
        xlabel('GDP per capita (\$, thousand)');
        vec_pos = get(get(gca,'xlabel'),'Position');
        set(get(gca,'xlabel'),'Position',vec_pos + [0 -0.5 0]);
    end
    if (j==1 && k==2);
        ylabel('Societal Cost (\$, trillion)');
        vec_pos = get(get(gca,'ylabel'),'Position');
        set(get(gca,'ylabel'),'Position',vec_pos + [-0.15 -ylm/2 0]);
    end
    set(gca,'FontSize',fs);   
                         
    end
end

%% STRATEGY COMPARISON

sumstatt = readtable('sumstatt_pc.csv');
sumstatt = table2array(sumstatt);

igroups    = {'LLMIC','UMIC','HIC'};
diseases   = {'Influenza 2009','Influenza 1918',...
              'Covid Omicron','Covid Wildtype',...
              'SARS'};
strategies = {'No Closures','Economic Closures','Elimination'};     
optimal    = cell(length(diseases)*length(igroups),4+2*length(strategies));          
          
for j = 1:length(diseases);
    for i = 1:length(igroups);
    
    inp1   = igroups(i);
    inp2   = diseases(j);
    [~,ko] = min(sumstatt((j-1)*length(strategies)+1:j*length(strategies),i));
    inp3   = strategies(ko);
    
    Topt = readtable(strcat(string(inp1),'_',string(inp2),'_',string(inp3),'.csv'));
    gopt = 365*sum(table2array(Topt(:,331:375)),2); 
    yopt = 100*Topt.SEC./gopt;
    
    optimal{i+(j-1)*length(igroups),1} = inp2;
    optimal{i+(j-1)*length(igroups),2} = inp1;
    optimal{i+(j-1)*length(igroups),3} = inp3;    
    optimal{i+(j-1)*length(igroups),4} = mean(yopt(~isnan(yopt)));

    for k = 1:length(strategies);    
    Talt = readtable(strcat(string(inp1),'_',string(inp2),'_',string(strategies(k)),'.csv'));
    galt = 365*sum(table2array(Talt(:,331:375)),2); 
    yalt = 100*Talt.SEC./galt;

    ydif = yalt-yopt;
    cdif = mean(ydif(~isnan(ydif)));
    pval = sum(ydif<0)/numel(ydif);
    
    if k==ko;
        cdif = NaN;
        pval = NaN;
    end
    
    optimal{i+(j-1)*length(igroups),4+2*k-1} = cdif;
    optimal{i+(j-1)*length(igroups),4+2*k}   = pval;
    end
    
    end
end

optimal = cell2table(optimal);
writetable(optimal,'sumstatt_optimal.csv');

%% VIOLINS 2

f  = figure('Units','centimeters','Position',[0 0 26 15]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

diseases   = {'Influenza 2009','Influenza 1918',...
              'Covid Omicron','Covid Wildtype',...
              'SARS'};
strategies = {'No Closures','Economic Closures','Elimination'};
ylimmax    = [250,1200,350,450,600];

for j = 1:length(diseases);
    for k = 1:length(strategies);
    
    inp2 = diseases(j);
    inp3 = strategies(k);
    T1   = readtable(strcat('LLMIC_',string(inp2),'_',string(inp3),'.csv'));
    T2   = readtable(strcat('UMIC_',string(inp2),'_',string(inp3),'.csv'));
    T3   = readtable(strcat('HIC_',string(inp2),'_',string(inp3),'.csv'));

    gdp1 = 365*sum(table2array(T1(:,331:375)),2);
    gdp2 = 365*sum(table2array(T2(:,331:375)),2);
    gdp3 = 365*sum(table2array(T3(:,331:375)),2);    
    y1   = 100*T1.SEC./gdp1;
    y2   = 100*T2.SEC./gdp2;
    y3   = 100*T3.SEC./gdp3;
    
    yran1 = linspace(0,max(y1),1000);
    yran2 = linspace(0,max(y2),1000);
    yran3 = linspace(0,max(y3),1000);
    pdf1  = ksdensity(y1,yran1,'Support','positive');%,'Bandwidth',0.1);
    pdf2  = ksdensity(y2,yran2,'Support','positive');%,'Bandwidth',0.1);
    pdf3  = ksdensity(y3,yran3,'Support','positive');%,'Bandwidth',0.1);
    
    ymax = max([pdf1,pdf2,pdf3]);
    pdf1 = pdf1/ymax/2;
    pdf2 = pdf2/ymax/2;
    pdf3 = pdf3/ymax/2;
    
    pdf1(pdf1<0.001) = 0;
    pdf2(pdf2<0.001) = 0;
    pdf3(pdf3<0.001) = 0;
    pdp1             = pdf1;
    pdp2             = pdf2;
    pdp3             = pdf3;
    pdp1(pdp1==0&circshift(pdp1,1)==0&circshift(pdp1,-1)==0) = NaN;
    pdp2(pdp2==0&circshift(pdp2,1)==0&circshift(pdp2,-1)==0) = NaN;
    pdp3(pdp3==0&circshift(pdp3,1)==0&circshift(pdp3,-1)==0) = NaN;
    
    lp1 = mean(T1.VLYL(~isnan(T1.SEC))./T1.SEC(~isnan(T1.SEC)));
    lp2 = mean(T2.VLYL(~isnan(T2.SEC))./T2.SEC(~isnan(T2.SEC)));
    lp3 = mean(T3.VLYL(~isnan(T3.SEC))./T3.SEC(~isnan(T3.SEC)));
    gp1 = mean(T1.GDPL(~isnan(T1.SEC))./T1.SEC(~isnan(T1.SEC)));
    gp2 = mean(T2.GDPL(~isnan(T2.SEC))./T2.SEC(~isnan(T2.SEC)));
    gp3 = mean(T3.GDPL(~isnan(T3.SEC))./T3.SEC(~isnan(T3.SEC)));
       
    s = subplot(length(strategies),length(diseases),j+(k-1)*length(diseases));
    d = 0.2;
    i = 0.0125+d*(0.975*1/5)+(j-1)*(0.975*1/5);
    w = (1-d)*(0.975*1/5);
    u = (d/2)*(1/3)+(3-k)*(1/3);
    h = (1-d)*(1/3);
    set(s,'position',[i u w h]);
    hold on;
    
    col1 = [1 .5 0];
    col2 = [0 0 1];
    col3 = [1 1 0];
    
    xval = 0.5;
    yran = yran1;
    ydat = y1;
    pdf  = pdf1;
    pdp  = pdp1;
    lp   = lp1;
    gp   = gp1;
    plot(xval+lp*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-lp*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval-fliplr(lp*pdf)],[yran,fliplr(yran)],col1,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval+fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-lp*pdf,xval-fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+(lp+gp)*pdf,xval+fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-(lp+gp)*pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    xval = 1.5;
    yran = yran2;
    ydat = y2;
    pdf  = pdf2;
    pdp  = pdp2;
    lp   = lp2;
    gp   = gp2;
    plot(xval+lp*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-lp*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval-fliplr(lp*pdf)],[yran,fliplr(yran)],col1,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval+fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-lp*pdf,xval-fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+(lp+gp)*pdf,xval+fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-(lp+gp)*pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    xval = 2.5;
    yran = yran3;
    ydat = y3;
    pdf  = pdf3;
    pdp  = pdp3;
    lp   = lp3;
    gp   = gp3;
    plot(xval+lp*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-lp*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval-fliplr(lp*pdf)],[yran,fliplr(yran)],col1,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-(lp+gp)*pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+lp*pdf,xval+fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-lp*pdf,xval-fliplr((lp+gp)*pdf)],[yran,fliplr(yran)],col2,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+(lp+gp)*pdf,xval+fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    fill([xval-(lp+gp)*pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col3,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    
    xlim([0 3]);
    ylm = ylimmax(j);
    ylim([0 ylm]);
    xticks([0.5:1:2.5]);
    yticks([0 ylm/2 ylm]);
    grid on;
    grid minor;
    box on;
    xticklabels({''});
    if (j==1 && k==2);
        ylabel('Societal Cost (\% of GDP)');
        vec_pos = get(get(gca,'ylabel'),'Position');
        set(get(gca,'ylabel'),'Position',vec_pos + [-0.125 0 0]);
    end
    set(gca,'FontSize',fs);   
                         
    end
end

%% TORNADO PLOTS

f  = figure('Units','centimeters','Position',[0 0 26 15]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

diseases   = {'Influenza 2009','Influenza 1918',...
              'Covid Omicron','Covid Wildtype',...
              'SARS'};
strategies = {'No Closures','Economic Closures','Elimination'};

infoplot = readtable('infoplot.csv');
infoplot = table2array(infoplot);
infoplot(:,1)         = [];
infoplot(:,2)         = [];
infoplot(5*8+1:6*8,:) = [];
infoplot(1*8+1:2*8,:) = [];

comboso  = categorical({'Demography','Surveillance','Distancing','Bed Capacity','Vaccination'});
comboso  = reordercats(comboso,cellstr(comboso)');
colourso = [[0 0 0];[0 1 0];[1 0 0];[.58 0 .83];[1 .74 .53]];

for j = 1:length(diseases);
    for k = 1:length(strategies);

    infos      = infoplot((8*(j-1)+1):(j*8),k);
    infos(2:4) = [];
    [~,inds]   = sort(infos,'ascend');
    combos     = reordercats(comboso,inds);
    colours    = colourso(inds,:);
    
    s = subplot(length(strategies),length(diseases),j+(k-1)*length(diseases));
    d = 0.2;
    i = d*(0.975*1/5)+(j-1)*(0.975*1/5);
    w = (1-d)*(0.975*1/5);
    u = 0.005+(d/2)*(1/3)+(3-k)*(1/3);
    h = (1-d)*(1/3);
    set(s,'position',[i u w h]);
    hold on;
    
    b           = barh(combos,infos);
    b.FaceColor = 'flat';
    b.CData     = colours;  
    
    xlim([0 0.70]);
    %ylim([]);
    xticks([0:0.35:0.70]);
    grid on;
    grid minor;
    box on;
    %xticklabels({});
    yticklabels({''});
    set(gca,'FontSize',fs);    
%     if j==5&&k==1;
%     hold on;
%     hBB=barh(combos,nan(5,5));                % now create the dummy four bar handles
%     hLG=legend(hBB,(combos),'location','southeast');  
%     for k=1:5;
%         hBB(k).FaceColor     = colourso(k,:);
%     end
%     end

    end
end

%%

filename = 'VOI_Covid Wildtype_Economic Closures.csv';
T        = readtable(filename);
gdp      = 365*sum(table2array(T(:,331:375)),2);       
SEC      = 100*T.SEC./gdp;
VLYL     = 100*T.VLYL./gdp;
VSYL     = 100*T.VSYL./gdp;
GDPL     = 100*T.GDPL./gdp;

f      = figure('Units','centimeters','Position',[0 0 40 20]);
params = {'Tres','t_tit','trate','sdl','sdb','Hmax','t_vax','arate','puptake'};
output = SEC;
for i = 1:length(params);
    
    param = params{i};
    %if strcmp(param,'sdb');
    %    param = log(T.(param));
    %else
        param = T.(param);
    %end    
    R     = corrcoef(param,output,'rows','complete');
    subplot(ceil(sqrt(length(params))),ceil(sqrt(length(params))),i);
    scatter(param,output,'.');
    title(round(R(1,2),2));
    xlabel(params{i});
    
end

%% P2

f  = figure('Units','centimeters','Position',[0 0 26 6]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
fs = 12;

igroups    = {'LLMIC'};
diseases   = {'Influenza 2009','Influenza 1918',...
              'Covid Omicron','Covid Wildtype',...
              'SARS'};
strategies = {'No Closures'};
%ylimmax    = [5000,150,10,20,250];
%ylimmin    = [-10,0,0,0,0];
%ylimmax    = [10 50 50 50 25];

for j = 1:4;%length(diseases);
    
    inp1 = igroups(1);
    inp2 = diseases(j);
    inp3 = strategies(1);
    T    = readtable(strcat('P2L',string(inp1),'_',string(inp2),'_',string(inp3),'.csv'));

    gdp  = 365*sum(table2array(T(:,331:375)),2);  
%     y0   = 100*T.SEC./gdp;
%     yt   = 100*T.SECt./gdp;
    y0   = T.SEC;
    yt   = T.SECt;
    sben = y0-yt;%societal benefit in $mill    
    intcost= T.Cost;
    %nben = sben-p2cost;

    ydat = nben;
    yran = linspace(min(ydat),max(ydat),1000);
    %pdf  = ksdensity(ydat(ydat>0),yran,'Support','positive');%,'Bandwidth',0.1);
    pdf  = ksdensity(ydat,yran);%,'Bandwidth',0.1);
    
    ymax = max(pdf);
    pdf  = pdf/ymax/2;
    
    pdf(pdf<0.001) = 0;
    pdp            = pdf;
    pdp(pdp==0&circshift(pdp,1)==0&circshift(pdp,-1)==0) = NaN;
          
    s = subplot(length(strategies),length(diseases),j);
    d = 0.2;
    i = 0.0125+d*(0.975*1/5)+(j-1)*(0.975*1/5);
    w = (1-d)*(0.975*1/5);
    u = (d/2)*(1/1)+(1-1)*(1/1);
    h = (1-d)*(1/1);
    set(s,'position',[i u w h]);
    hold on;
    
    xval = 0.5;
    col  = [1,1,0];%[0 0 1];%[1 .5 0]
    plot(xval+pdp,yran,'Linewidth',0.5,'Color','black');
    plot(xval-pdp,yran,'Linewidth',0.5,'Color','black');
    fill([xval+pdf,xval-fliplr(pdf)],[yran,fliplr(yran)],col,'FaceAlpha',0.5,'EdgeColor','none');
    plot(xval,quantile(ydat,0.50),'Marker','s','MarkerSize',7.5,'Linewidth',1.5,'MarkerEdgeColor','black','MarkerFaceColor','white');
    
    xlim([0 1]);
    ylim([min(ydat) max(ydat)]);
    xticks([0.5]);
    %yticks([0 ylm/2 ylm]);
    grid on;
    grid minor;
    box on;
    xticklabels({''});
    if (j==1);
        %ylabel('Societal Cost (\$, trillion)');
        %ylabel('Societal Benefit (\% of GDP)');
        %vec_pos = get(get(gca,'ylabel'),'Position');
        %set(get(gca,'ylabel'),'Position',vec_pos + [-0.05 0 0]);
    end
    set(gca,'FontSize',fs);  
    
end
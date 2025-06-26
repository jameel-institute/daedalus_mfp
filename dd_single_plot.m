function f=dd_single_plot(data,g1,p2,f1,cost,inp1,inp2,inp3)

if strcmp(inp3, 'School Closures')
    inp3 = 'Reactive/Sustained-School Closures';
elseif strcmp(inp3, 'Economic Closures')
    inp3 = 'Reactive/Reactive-School Closures';
end

%for schematic figure s4
% add to last no closure event tflag: + max(0,(3*365)-t);
% add legend 2 to top subfigure only
% remove dashes in losses barchart for subfig 2 and fix ylmt
% manually specify income group of country in legend 1
% dd_single('Rwanda','Influenza 1957','No Closures');
% dd_single('Argentina', 'Influenza 1918', 'School Closures');
% dd_single('Thailand', 'Covid Omicron', 'Economic Closures');
% dd_single('United Kingdom', 'Covid Wildtype', 'Elimination');

ln        = length(data.NNs);
lx        = length(data.obj);
tvec      = data.tvec;
thresh    = p2.Hmax;

t1    = g1(:,1); 
I1    = g1(:,2);
h1    = g1(:,3);
d1    = sum(g1(:,16:19),2);
% ddiff = diff(d1,1);
% tdiff = diff(t1,1);
% d1    = [0;ddiff./tdiff];
asc_a = g1(:,5);
asc_s = g1(:,6);
beta  = g1(:,7);
v1    = g1(:,8);
v2    = g1(:,9);
v3    = g1(:,10);
v4    = g1(:,11);

isoasyu = f1(:,1+2*lx+0*ln+[1:ln]);
isoasyv = f1(:,1+2*lx+1*ln+[1:ln]);
isosymu = f1(:,1+2*lx+2*ln+[1:ln]);
isosymv = f1(:,1+2*lx+3*ln+[1:ln]);
Q       = sum(isoasyu+isoasyv+isosymu+isosymv,2);

scal4 = sum(data.Npop)/(10^4);
scal5 = sum(data.Npop)/(10^5);
scal6 = sum(data.Npop)/(10^6);
scal7 = sum(data.Npop)/(10^7);
scal8 = sum(data.Npop)/(10^8);
scal9 = sum(data.Npop)/(10^9);
%maxY  = max([100000,(5/4)*d1'/scal5,(5/4)*thresh/scal4,(5/4)*h1'/scal4,(5/4)*I1'/scal3]);
maxY  = 200000;%ceil(max([d1'/scal5,thresh/scal4,h1'/scal4,I1'/scal3])/10000)*10000;

T               = repmat(t1',lx+1,1);
S               = 0.5:1:lx+0.5;
x               = f1(:,1+[1:lx])';
X               = [x;ones(1,length(t1))];
notEd           = [1:(data.EdInd-1),(data.EdInd+1):lx];

%% EPIDEMIC TRAJECTORY

f  = figure('Units','centimeters','Position',[0 0 21 9]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultAxesColorOrder',[[0 0 0];[0 0 0]]);
set(f,'DefaultAxesFontSize',12);
fs = 12;
lw = 2;

ax = gca;
ax.Position = [0.05 0.20 0.90 0.70];
hold on;

for i = 1:length(tvec)-1; 
    a     = [tvec(i) tvec(i) tvec(i+1) tvec(i+1)];
    b     = [0 maxY maxY 0];
    pind  = tvec(i)<t1 & t1<tvec(i+1);
    px    = x(notEd,pind);%this is all sectors excl education
    ps    = x(data.EdInd,pind);
    alpha = 1-mean(px.^3,'all');
    if alpha == 0;
        continue;
    end
    fill(a,b,'yellow','facealpha',alpha,'linewidth',0.01,'EdgeColor','k');

    if mean(ps) < 1;
    
    x_width = 500;                        
    slope   = maxY / x_width;           
    spacing = 50;            
    x_start = tvec(i);
    x_end   = tvec(i+1);
    y_bot   = 0;
    y_top   = maxY;
    for x0 = x_start - y_top / slope : spacing : x_end
        x1 = x0 + y_top / slope;
        x_vals = [x0, x1];
        y_vals = [y_bot, y_top];
        if x_vals(2) < x_start || x_vals(1) > x_end
            continue;
        end
        x_vals = max(min(x_vals, x_end), x_start);
        y_vals = y_bot + (x_vals - x0) * slope;
        plot(x_vals, y_vals, 'b', 'LineWidth', 0.5);
    end
    end

end
%plot([p2.Tres,p2.Tres],[0,maxY],'-','linewidth',0.01,'color',0.5*[1,1,1]);

yyaxis left;
hh3 = plot([0,tvec(end)],[thresh,thresh]/scal8,'--','linewidth',lw,'color',[0.7, 0, 0.9]);
hh4 = plot(t1,d1/scal9,'-','linewidth',lw,'color','red');
hh2 = plot(t1,h1/scal8,'-','linewidth',lw,'color',[0.7, 0, 0.9]);
%hh1 = plot(t1,Q/scal7,'-','linewidth',lw,'color','green');
hh1 = plot(t1,(I1)/scal7,':','linewidth',lw,'color','black');
%hh0 = plot(t1,beta*maxY,'-','linewidth',lw,'color','blue');
%hh0 = plot(t1,asc_s/scal8,'-','linewidth',lw,'color','green');
yyaxis right;
hh5 = plot(t1,100*(v1+v2+v3+v4)/sum(data.Npop),'-','linewidth',lw,'color',[1, 0.6, 0.4]);
xlim([1 3*365]);
set(gca,'xtick',[[1,91,182,274],...
             365+[1,91,182,274],...
           2*365+[1,91,182,274],...
           3*365+[1,91,182,274]]);
set(gca,'xticklabels',{'Year 1','Apr','Jul','Oct',...
                       'Year 2','Apr','Jul','Oct',...
                       'Year 3','Apr','Jul','Oct',...
                       'Year 4','Apr','Jul','Oct'});
xtickangle(45);
yyaxis left;
ylim([0 maxY]);
ax = gca;
ax.YAxis(1).Exponent = 3;
ax.YColor = 'k';
yyaxis right;
ylim([0 100]);
ylabel('\%');
set(get(gca,'ylabel'),'Rotation',0);
vec_pos = get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos + [-60 60 0]);
ax = gca;
ax.YColor = 'k';

grid on;
box on;
% legend([hh1,hh2,hh3,hh4,hh5],'Prevalence (per 10m)','Hospital Occupancy (per 100m)','Hospital Capacity (per 100m)',... %%%%%
%                          'Daily Deaths (per 1b)','Vaccine Coverage (\%)','location','northeast');
t = text(8,86.6,strvcat(inp1,strrep([inp2 ' X'], ' ', '-'),inp3),'FontSize',fs);
t.BackgroundColor = [1 1 1];
t.EdgeColor       = [0 0 0];
set(gca,'FontSize',fs);
yyaxis left;

%% MITIGATION MEASURES

% f = figure('Units','centimeters','Position',[0 0 10 10]);
% set(f,'defaulttextInterpreter','latex');
% set(f,'defaultAxesTickLabelInterpreter','latex');
% set(f,'defaultLegendInterpreter','latex');
% set(f,'defaultColorbarTickLabelInterpreter','latex');
% set(f,'DefaultAxesFontSize',12);
% 
% ax = gca;
% ax.Position = [0.175 0.20 0.80 0.70];
% 
% h = pcolor(T,S,X);
% 
% axis([-30,3*365,0.5,lx+0.5]);
% %axis square;
% set(gca,'layer','top');
% set(gca,'xtick',[[1,182],...
%              365+[1,182],...
%            2*365+[1,182]]);
% set(gca,'ytick',0.5+[5:5:lx]);
% grid on;
% box on;
% set(gca,'xticklabels',{'Year 1 - Jan','Jul',...
%                        'Year 2 - Jan','Jul',...
%                        'Year 3 - Jan','Jul'});       
% set(gca,'yticklabels',{'5','10','15','20','25','30','35','40','45'});
% xtickangle(45);
% %xlabel('Time');
% ylabel('Sector')
% set(h,'EdgeColor','none');
% set(gca,'FontSize',fs);
% 
% colormap(hot);
% caxis([0,1]);
% colorbar;

%% COSTS

f = figure('Units','centimeters','Position',[0 0 6 9]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);

ax          = gca;
ax.Position = [0.05 0.20 0.65 0.70];

labs = categorical(["VLYL","GDPL","VSYL"]);
labs = reordercats(labs,[2 1 3]);
y    = 100*[cost(3,lx+1)   cost(3,lx+2)   sum(cost(3,[1:lx,lx+3])) cost(3,ln);...
            sum(cost(4,:)) sum(cost(5,:)) 0                        0;...
            cost(8,lx+1:ln)]/sum(365*data.obj);
b    = bar(labs,y,'stacked','FaceColor','flat');
for i = 1:length(b)
    b(i).EdgeColor = 'none';
end

xtickangle(45);
ymax = 100*max([sum(cost(3,:)),sum(cost(4:5,:),'all'),sum(cost(8,:))])/sum(365*data.obj);
ylmt = ceil(ymax/10)*10;
ylmt = ylmt + 10*((ylmt-ceil(ymax))<(0.25*ylmt));%%%%%
ylim([0 ylmt]);
ylh = ylabel('Socioeconomic Loss (\% of GDP)');
ax = gca;
ax.YAxisLocation = 'right';
ylh.Rotation = 270;
pos = get(ylh, 'Position');
set(ylh, 'Position', pos + [0.85 0 0]);  
%text(0.6,ylmt*1.05,'\%','FontSize',fs);

grid on;
grid minor;
box on;
b(1).CData = [1.00 0.00 0.00;...   
              1.00 1.00 0.00;...   
              0.00 0.00 1.00];     
b(2).CData = [1.00 0.00 0.00;...
              1.00 1.00 0.00;...
              0.00 0.00 1.00];
b(3).CData = [1.00 0.00 0.00;...
              1.00 0.84 0.00;...
              0.00 0.00 1.00];
b(4).CData = [1.00 0.00 0.00;...
              1.00 1.00 0.00;...
              0.00 0.00 1.00];
set(gca,'FontSize',fs);

ax.YTick(end) = [];%%%%%
axes('Position',[.665 .81 .05 .05]);
px     = [1 5];
py1    = [1 2];
height = 1;
py2    = py1 + height;
plot(px,py1,'k','LineWidth',2);hold all;
plot(px,py2,'k','LineWidth',2);hold all;
fill([px flip(px)],[py1 flip(py2)],'w','EdgeColor','none');
box off;
axis off;

end
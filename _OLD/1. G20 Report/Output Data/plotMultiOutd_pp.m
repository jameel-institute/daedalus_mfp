function f=plotMultiOutd_pp(f1,x1,tvec,data,pr)

numSectors=length(data.G);
numPeriods=length(tvec)-1;
dodiff=1;

%%

f=figure('Units','centimeters','Position',[0 0 30 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'DefaultAxesFontSize',12);fs=12; lw=2;
ax=gca;
ax.Position=[0.05 0.175 0.90 0.75];
cmap=lines(2);
thresh=pr.Hmax;
numThresh=length(thresh);
t1=f1(:,1); 
s1=f1(:,2); 
I1=f1(:,3);
h1=f1(:,4);
d1=f1(:,5);
v1=f1(:,6);%heSimCovid19 output
v2=f1(:,7);%heSimCovid19 output
v3=f1(:,8);%heSimCovid19 output
v4=f1(:,9);%heSimCovid19 output
if dodiff==1
    inc1=-diff(s1,1);%f1=Sout
    tdiff=diff(t1,1);
    inc1=inc1./tdiff;%repmat - older version?
end
hold on

scal1=sum(data.Npop)/(10^5);
scal2=sum(data.Npop)/(10^6);
scal3=sum(data.Npop)/(10^7);
scal4=sum(data.Npop)/(10^8);

%maxY=max([5*max(I1/scal1)/4,5*max(d1)/4,100000,pr.Hmax,5*max(h1)/4]);%inc2;h2
maxY=max([100000,(5/4)*d1'/scal4,(5/4)*thresh/scal4,(5/4)*h1'/scal4,(5/4)*I1'/scal3]);%inc2;h2

plot(tvec(2)*[1,1],[0,maxY],'k-','linewidth',0.01);

for i=3:2:length(tvec)-1 
    %for i=2:length(tvec)-1
    %plot(tvec(i)*[1,1],[0,maxY],'k--','linewidth',1)
    a=[tvec(i) tvec(i) tvec(i+1) tvec(i+1)];
    b=[-100 maxY+100 maxY+100 -100];
    fill(a,b,'cyan','linewidth',0.01,'facealpha',.2);
end

for j=1:numThresh
    plot([0,tvec(end)],[thresh(j),thresh(j)]/scal4,':','linewidth',lw,'color',.5*[1,1,1])
end

%hh1=plot(t1(1:end-1)+0.5*tdiff,inc1,'-','linewidth',lw,'color','yellow');
hh5=plot(t1,(v1+v2+v3+v4)/scal1,'-','linewidth',lw,'color','green');
hh4=plot(t1,d1/scal4,'-','linewidth',lw,'color','black');
hh3=plot(t1,h1/scal4,'-','linewidth',lw,'color','magenta');
hh2=plot(t1,I1/scal3,'-','linewidth',lw,'color','red');

points=tvec+5;
pointsy=.95*maxY;
txt={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30',...
     '31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50','51','52','53','54','55','56','57','58','59','60'};
text(5,pointsy,'PRE','fontsize',fs);
%text(tvec(2)+5,pointsy,'LD','fontsize',15);
for i=2:numPeriods
    text(points(i),pointsy,txt{i-1},'fontsize',fs);
end
axis ([1,365*3,0,maxY])
xlabel('Time','FontSize',fs);
ylabel('Number','FontSize',fs);%yvar
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-10 0 0]);
set(gca,'FontSize',fs);

%xticks([1,32,61,92,122,153,183,214,245,275,306,336,367,398])
%xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb'})
set(gca,'xtick',[[1,32,60,91,121,152,182,213,244,274,305,335],365+[1,32,60,91,121,152,182,213,244,274,305,335],2*365+[1,32,60,91,121,152,182,213,244,274,305,335]]);
set(gca,'xticklabels',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',...
                       'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec',...
                       'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
xtickangle(45);
ax = gca;
ax.YAxis.Exponent = 3;
%legend([hh1,hh2],'Inc.','Hosp. occ.','location','west')
legend([hh2,hh3,hh4,hh5],'Prevalence (per 10m)','Hospital Occupancy (per 100m)','Deaths (per 100m)','Vaccinated (per 100k)','location','northeast');%'Position',[-0.29 0.27 1 1]);
%legend([hh1,hh2,hh3,hh4],'Inc. (xmin)','HC (xmin)','Inc. (open)','HC (open)','location','west')
grid on;
%grid minor
box on
hold off

%%

f=figure('Units','centimeters','Position',[0 0 10 10]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax=gca;
ax.Position=[0.15 0.10 0.80 0.80];

if pr.sw==0
    x1=[ones(numSectors,1);x1];
    x=reshape(x1,numSectors,numPeriods);
else
    x1=repmat(reshape(x1,numSectors,2),1,length(tvec));
    x=[ones(numSectors,1),x1(:,1:length(tvec)-2)];
end

xvec=(1:numPeriods)';xvec=xvec-.5;
xlabels2=({'PRE','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20'});
yvec=(1:5:numSectors)';
ylab=num2str(yvec);
yvec=yvec-.5;
colormap hot;
imagesc(x)
set(gca,'YDir','normal')
xlabel('Period')
ylabel('Sector')
vec_pos=get(get(gca,'ylabel'),'Position');
set(get(gca,'ylabel'),'Position',vec_pos+[-0.5 0 0]);
set(gca,'fontsize',fs,'xtick',xvec,'xticklabels',xlabels2,'ytick',yvec,'yticklabels',ylab);%{'PRE','LD',xlab(3:end,:)}
axis square;%([0,numPeriods,.5,63.5])%yvec(1),yvec(end)+1])
xtickangle(45);
caxis([0,1])
colorbar
grid on
box on

end
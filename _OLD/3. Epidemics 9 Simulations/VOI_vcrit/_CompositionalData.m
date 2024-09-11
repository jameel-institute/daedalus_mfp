clist       = dir(fullfile('*.mat'));
clist       = {clist.name};
clist(end)  = [];
X           = zeros(length(clist),21);
for i = 1:length(clist);
    load(clist{i},'data');  
    X(i,:)  = data.Npop'/sum(data.Npop);  
end
X(:,19)     = sum(X(:,19:end),2);%zero values cause problems for fitting
X(:,20:end) = [];
writematrix(X,'Npop.csv');

d1 = sum(X(:,1:4),2);
%d2 = sum(X(:,5:13),2);
d3 = sum(X(:,14:end),2);
dx = d3 + d1/2;
dy = sqrt(3)*d1/2;

x1       = [0:0.01:1];
%x2       = [0:0.01:1];
x3       = [0:0.01:1];
[X1,X3]  = meshgrid(x1,x3);
X2       = 1 - X1 - X3;
X1       = X1(:);
X2       = X2(:);
X3       = X3(:);
inds     = find(X1+X3>1);
X1(inds) = [];
X2(inds) = [];
X3(inds) = [];
XX       = X3 + X1/2;
YY       = sqrt(3)*X1/2;
tri      = delaunay(XX,YY);

alp   = readmatrix('alpha.csv')';%[6 12 4];
alp   = [sum(alp(1:4)),sum(alp(5:13)),sum(alp(14:end))];%convolution of dirichlet is dirichlet with summed parameters
bnorm = exp(sum(gammaln(alp))-gammaln(sum(alp)));
F     = real((X1.^(alp(1)-1).*X2.^(alp(2)-1).*X3.^(alp(3)-1))/bnorm);
%fun   = @(x1,x2) (x1+x2<=1).*real((x1.^(alp(1)-1).*x2.^(alp(2)-1).*(1-x1-x2).^(alp(3)-1))/bnorm);
%test  = integral2(fun,0,1,0,1);
%fun   = betapdf(x1,alp(1),sum(alp)-alp(1));%marginals of dirichlet are beta distributed
%test  = trapz(x1,fun);

f = figure('Units','centimeters','Position',[0 0 9 9]);
set(f,'defaulttextInterpreter','latex');
set(f,'defaultAxesTickLabelInterpreter','latex');
set(f,'defaultLegendInterpreter','latex');
set(f,'defaultColorbarTickLabelInterpreter','latex');
set(f,'DefaultAxesFontSize',12);
ax = gca;
ax.Position = [0.05 0.15 0.70 0.70];
hold on;

trisurf(tri,XX,YY,F);
plot3(dx,dy,10^6*ones(size(dy)),'r.');
t  = [0:0.1:0.9];
l  = 0.01;
lx = l*cos(pi/3);
ly = l*sin(pi/3);
plot3([-l+t/2;1-t/2],    [sqrt(3)*t/2;sqrt(3)*t/2],          10^6*ones(2,length(t)),':','color',[0.5 0.5 0.5]);
plot3([lx+(1-t);(1-t)/2],[-ly*ones(size(t));sqrt(3)*(1-t)/2],10^6*ones(2,length(t)),':','color',[0.5 0.5 0.5]);
plot3([lx+(1+t)/2;t],    [ly+sqrt(3)*(1-t)/2;zeros(size(t))],10^6*ones(2,length(t)),':','color',[0.5 0.5 0.5]);
plot3([0,0.5,1;0.5,1,0], [0,sqrt(3)/2,0;sqrt(3)/2,0,0],      10^6*ones(2,3),            'color','black');

axis off;
for i = 1:length(t);
    ht = text(-0.07+t(i)/2,sqrt(3)*t(i)/2,num2str(100*t(i)));
    set(ht,'Rotation',0);
    ht = text(1-t(i),-0.02,num2str(100*t(i)));
    set(ht,'Rotation',-60);
    ht = text(0.01+(1+t(i))/2,0.01+sqrt(3)*(1-t(i))/2,num2str(100*t(i)));
    set(ht,'Rotation',60);
end
ht  = text(0.03,0.36,'School-Age (\%)');
set(ht,'Rotation',60);
ht  = text(0.3,-0.15,'Working-Age (\%)');
set(ht,'Rotation',0);
ht  = text(0.77,0.68,'Retired-Age (\%)');
set(ht,'Rotation',-60);
axis image;
view([0 90]);
colormap(flipud(bone));
shading interp;
cb  = colorbar;
pos = cb.Position;
set(cb,'Position',pos + [0.25 -0.07 0 0.17]);

rv        = gamrnd(alp,1);
rv        = rv./sum(rv);%see properties of dirichlet distribution
data.Npop = sum(data.Npop)*rv';

ads      = sum(nn(5:13));
wf       = data.NNs(1:45);
cap      = min(ads/sum(wf),1);
wf       = cap*wf;
data.NNs = [wf;nn(1);sum(nn(2:4));ads-sum(wf);sum(nn(14:end))];
%the assumpion is that the number in the workforce is the same, regardless of the adult population
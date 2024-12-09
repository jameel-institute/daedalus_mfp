function [f,data]=p2MakeDs(data,NN,x,hw)

%define number of strata, sectors, adult index, population by age (compressed to 16)
ln       = length(NN);
lx       = length(x);
adInd    = 3;
Npop     = data.Npop;
Npop(16) = sum(Npop(16:end));
Npop     = Npop(1:16);

%initialise contact matrices: notation consistent with SI of Haw et al. (2022)
matAL = zeros(ln,ln);%household
matAS = zeros(ln,ln);%school
matAH = zeros(ln,ln);%hospitality
matAT = zeros(ln,ln);%transport
matB  = zeros(ln,ln);%worker-worker
matC  = zeros(ln,ln);%consumer-worker

%% matAL: household

%start with 16x16 contact matrix, compress it to 4x4 and normalise it
CM     = data.CM;
CM     = [CM(:,1),sum(CM(:,2:4),2),sum(CM(:,5:13),2),sum(CM(:,14:16),2)];%sum of the columns
CM     = [CM(1,:);
          Npop(2:4)'*CM(2:4,:)/sum(Npop(2:4));
          Npop(5:13)'*CM(5:13,:)/sum(Npop(5:13));
          Npop(14:16)'*CM(14:16,:)/sum(Npop(14:16))];%weighted average of the rows
CMnorm = ([Npop(1),sum(Npop(2:4)),sum(Npop(5:13)),sum(Npop(14:16))]/sum(Npop))*sum(CM,2);%norm of the matrix
CM     = CM/CMnorm;

%expand CM into matAL, broadcast contacts from adults (columns) proportional to adult population
adRow                    = CM(adInd,:);
matAL(lx+1:end,lx+1:end) = CM;
matAL(1:lx,lx+1:end)     = repmat(adRow,lx,1);
NNrel                    = NN([1:lx,lx+adInd])/sum(NN([1:lx,lx+adInd]));%adult population proportion vector
matAL(:,[1:lx,lx+adInd]) = repmat(matAL(:,lx+adInd),1,lx+1).*repmat(NNrel',ln,1);

%specify magnitude of household contacts
ALmag = data.comm;
matAL = ALmag*matAL;

%% matAS: school

%school contacts are quadratic in x since we don't move people between strata
matAS(lx+1,lx+1) = data.schoolA1*x(data.EdInd)^2;
matAS(lx+2,lx+2) = data.schoolA2*x(data.EdInd)^2;

%% matAH: hospitality

%hospitality contacts split in proportion to total population (incl. pre-school), quadratic in x/psub since we don't move people between strata
NNrep = NN'/sum(NN);%total population proportion vector
psub  = data.NNs(data.HospInd);
psub  = sum(psub.*x(data.HospInd))/sum(psub);%constant from 0-1, weighted measure of how much sectors are open

matAH(lx+2,:)            = data.hospA2*NNrep*psub^2;
matAH([1:lx,lx+adInd],:) = data.hospA3*repmat(NNrep,lx+1,1)*psub^2;
matAH(ln,:)              = data.hospA4*NNrep*psub^2;

%% matAT: transport

%transport contacts split in proportion to workforce population, linear in x but quadratic in wfh
NNrea = repmat(NN(1:lx)'/sum(NN(1:lx)),lx,1);%workforce population proportion matrix

matAT(1:lx,1:lx) = data.travelA3.*NNrea.*repmat(x',lx,1).*repmat(1-hw,lx,1).*repmat(1-hw',1,lx);

%% matB: worker-worker

%worker-worker contacts linear in x but quadratic in wfh
matB(1:lx,1:lx) = diag(data.B.*x'.*(1-hw).*(1-hw));

%% matC: consumer-worker

%consumer-worker contacts linear in x but quadratic in wfh
matC(1:lx,:) = data.C'.*repmat(NNrep,lx,1).*x.*(1-hw');

%% matB and matC: workplace

%resulting number of workplace contacts
if ~isfield(data,'wnorm');
    data.wnorm = dot(sum(matB+matC,2),NN)/sum(NN([1:lx,lx+adInd]));
end

%specify magnitude of workplace contacts
matB = (data.workp/data.wnorm)*matB;
matC = (data.workp/data.wnorm)*matC;

%%

f = matAL + matAS + matAH + matAT + matB + matC;

end
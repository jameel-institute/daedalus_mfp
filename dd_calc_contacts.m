function f = dd_calc_contacts(data,x,hw)

%define number of strata, sectors, adult index and population
ln    = length(data.NNs);
lx    = length(x);
adInd = 3;
NN    = data.NNs;
NNtot = NN'/sum(NN);%total population proportion row vector
NNad  = NN([1:lx,lx+adInd])'/sum(NN([1:lx,lx+adInd]));%adult population proportion row vector

%initialise contact matrices: notation consistent with SI of Haw et al. (2022)
matAL  = zeros(ln,ln);%household
matAHT = zeros(ln,ln);%other
matAS  = zeros(ln,ln);%school
matB   = zeros(ln,ln);%worker-worker
matC   = zeros(ln,ln);%consumer-worker

%% matAL: household

%read in 4x4 matrix
AL = data.matAL;

%expand AL into matAL, broadcast contacts from adults (columns) proportional to adult population
matAL(lx+1:end,lx+1:end) = AL;
adRow                    = AL(adInd,:);
matAL(1:lx,lx+1:end)     = repmat(adRow,lx,1);
matAL(:,[1:lx,lx+adInd]) = repmat(matAL(:,lx+adInd),1,lx+1).*repmat(NNad,ln,1);

%household contacts are independent of x and hw

%% matAHT: other

%read in 4x4 matrix
AHT = data.matAHT;

%expand AHT into matAHT, broadcast contacts from adults (columns) proportional to adult population
matAHT(lx+1:end,lx+1:end) = AHT;
adRow                     = AHT(adInd,:);
matAHT(1:lx,lx+1:end)     = repmat(adRow,lx,1);
matAHT(:,[1:lx,lx+adInd]) = repmat(matAHT(:,lx+adInd),1,lx+1).*repmat(NNad,ln,1);

%45% of other contacts assumed to be leisure, from https://journals.plos.org/plosmedicine/article/file?id=10.1371/journal.pmed.0050074&type=printable
%leisure contacts are quadratic in psub(x) 
psub   = dot(x(data.HospInd),data.NNs(data.HospInd))/sum(data.NNs(data.HospInd));%constant from 0-1, weighted measure of how much sectors are open
matAHT = 0.55*matAHT + 0.45*matAHT*psub^2;

%% matAS: school

%read in 4x4 matrix
AS = data.matAS;

%expand AS into matAS, broadcast contacts from adults (columns) proportional to adult population
matAS(lx+1:end,lx+1:end) = AS;
adRow                    = AS(adInd,:);
matAS(1:lx,lx+1:end)     = repmat(adRow,lx,1);
matAS(:,[1:lx,lx+adInd]) = repmat(matAS(:,lx+adInd),1,lx+1).*repmat(NNad,ln,1);

%school contacts are quadratic in x
matAS(lx+1:lx+2,lx+1:lx+2) = matAS(lx+1:lx+2,lx+1:lx+2)*x(data.EdInd)^2;

%% matB: worker-worker

%read in value and calculate unadjusted magnitude (in absence of NPIs)
workp = data.workp;%average number of workplace contacts per adult (working and non-working)
unmag = dot(data.B + data.C, NNad(1:lx));%consistent with definition above

%expand B into matB and scale
matB(1:lx,1:lx) = (workp/unmag)*diag(data.B);

%worker-worker contacts quadratic in x and hw
matB(1:lx,1:lx) = matB(1:lx,1:lx).*repmat(x',lx,1).*repmat(x,1,lx).*repmat(1-hw,lx,1).*repmat(1-hw',1,lx);

%% matC: consumer-worker

%expand C into matC and scale
matC(1:lx,:) = (workp/unmag)*data.C'.*repmat(NNtot,lx,1);

%consumer-worker contacts linear in x and hw
matC(1:lx,:) = matC(1:lx,:).*x.*(1-hw');

%%

f = matAL + matAHT + matAS + matB + matC;

end
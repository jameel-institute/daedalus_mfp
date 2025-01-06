function [Rt,ev] = dd_calc_Rt(dis,h,g2,S,Shv1,Sv1,N,D,betamod,sig1,sig2,sig3,sig4,tm_a,tm_s)

ln = length(N);
F  = zeros(10*ln,10*ln);

N(N==0) = 1;
ARssh   = dis.beta.*betamod.*D./repmat(N',ln,1).*repmat(S+Shv1,1,ln);
ARsv    = dis.beta.*betamod.*D./repmat(N',ln,1).*(1-dis.scv1).*repmat(Sv1,1,ln);
onesn   = ones(ln,1);

F(1:ln,     2*ln+1:end) = [dis.red.*ARssh,  ARssh,  tm_a.*dis.red.*ARssh,  tm_s.*ARssh, ...  
                           dis.red.*(1-dis.trv1).*ARssh,  (1-dis.trv1).*ARssh,  tm_a.*dis.red.*(1-dis.trv1).*ARssh,  tm_s.*(1-dis.trv1).*ARssh];
F(ln+1:2*ln,2*ln+1:end) = [dis.red.*ARsv,  ARsv,  tm_a.*dis.red.*ARsv,  tm_s.*ARsv, ...  
                           dis.red.*(1-dis.trv1).*ARsv,  (1-dis.trv1).*ARsv,  tm_a.*dis.red.*(1-dis.trv1).*ARsv,  tm_s.*(1-dis.trv1).*ARsv];

vvec                      = [(sig1+sig2+sig3+sig4).*onesn;
                             (sig1+sig2+sig3+sig4).*onesn;
                              dis.g1.*onesn;
                             (g2+h).*onesn;
                              dis.g1.*onesn;
                             (g2+h).*onesn;
                              dis.g1.*onesn;
                             (dis.g2_v1+dis.h_v1).*onesn;
                              dis.g1.*onesn;
                             (dis.g2_v1+dis.h_v1).*onesn];
V                         = diag(vvec);
V(2*ln+1:3*ln,1:ln)       = diag(-sig1.*onesn);
V(3*ln+1:4*ln,1:ln)       = diag(-sig2.*onesn);
V(4*ln+1:5*ln,1:ln)       = diag(-sig3.*onesn);
V(5*ln+1:6*ln,1:ln)       = diag(-sig4.*onesn);
V(6*ln+1:7*ln,ln+1:2*ln)  = diag(-sig1.*onesn);
V(7*ln+1:8*ln,ln+1:2*ln)  = diag(-sig2.*onesn);
V(8*ln+1:9*ln,ln+1:2*ln)  = diag(-sig3.*onesn);
V(9*ln+1:10*ln,ln+1:2*ln) = diag(-sig4.*onesn);

NGM = F/V;
if nargout == 1;
    Rt      = eigs(NGM,1,'largestreal');
else
    [ev,Rt] = eigs(NGM,1,'largestreal');
end
    
end
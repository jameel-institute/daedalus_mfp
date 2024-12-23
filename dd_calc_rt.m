function [rt,ev] = dd_calc_rt(dis,h,g2,S,Shv1,Sv1,N,D,betamod,sig1,sig2,sig3,sig4,tm_a,tm_s)
    
    N(N==0) = 1;
    ln      = length(N);
    F       = zeros(10*ln,10*ln);
    
    FOIu = dis.beta.*betamod.*D./repmat(N',ln,1).*repmat(S+Shv1,1,ln);
    FOIv = dis.beta.*betamod.*(1-dis.scv1).*D./repmat(N',ln,1).*repmat(Sv1,1,ln);
    
    F(1:ln,     2*ln+1:end) = [dis.red.*FOIu,  FOIu,  tm_a.*dis.red*FOIu,  tm_s.*FOIu, ...  
                               dis.red.*(1-dis.trv1).*FOIu,  (1-dis.trv1).*FOIu,  tm_a.*dis.red.*(1-dis.trv1).*FOIu,  tm_s.*(1-dis.trv1).*FOIu];
    F(ln+1:2*ln,2*ln+1:end) = [dis.red.*FOIv,  FOIv,  tm_a.*dis.red.*FOIv,  tm_s.*FOIv, ...  
                               dis.red.*(1-dis.trv1).*FOIv,  (1-dis.trv1).*FOIv,  tm_a.*dis.red.*(1-dis.trv1).*FOIv,  tm_s.*(1-dis.trv1).*FOIv];

    onesn = ones(ln,1);
    vvec  = [(sig1+sig2+sig3+sig4).*onesn;
             (sig1+sig2+sig3+sig4).*onesn;
              dis.g1.*onesn;
             (g2+h).*onesn;
              dis.g1.*onesn;
             (g2+h).*onesn;
              dis.g1.*onesn;
             (dis.g2_v1+dis.h_v1).*onesn;
              dis.g1.*onesn;
             (dis.g2_v1+dis.h_v1).*onesn];
    V     = diag(vvec);
    
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
        rt      = eigs(NGM,1,'largestreal');
    else
        [ev,rt] = eigs(NGM,1,'largestreal');
    end
    
end
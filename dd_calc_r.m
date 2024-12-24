function r = dd_calc_r(dis,h,g2,mu,g3,S,Shv1,Sv1,Ina,Ins,Isa,Iss,Inav1,Insv1,Isav1,Issv1,N,D,betamod,sig1,sig2,sig3,sig4,tm_a,tm_s)

    N(N==0) = 1;
    ln      = length(N);
    J       = zeros(18*ln,18*ln);
    
    I     = (dis.red*Ina+Ins) + tm_a*dis.red*(Isa+(1-dis.trv1)*Isav1) + (1-dis.trv1)*(dis.red*Inav1+Insv1) + tm_s.*(Iss+(1-dis.trv1)*Issv1);
    FOI   = dis.beta.*betamod.*D.*(repmat(I',ln,1)./repmat(N',ln,1));
    ARs   = dis.beta.*betamod.*D./repmat(N',ln,1).*repmat(S,1,ln);
    ARsh  = dis.beta.*betamod.*D./repmat(N',ln,1).*repmat(Shv1,1,ln);
    ARsv  = dis.beta.*betamod.*D./repmat(N',ln,1).*(1-dis.scv1).*repmat(Sv1,1,ln);
    onesn = ones(ln,1);
    
    %S
    J(1:ln,1:ln)                = -FOI;
    J(1:ln,2*ln+1:3*ln)         = diag(dis.nuv1.*onesn);
    J(1:ln,5*ln+1:13*ln)        = -[dis.red.*ARs,  ARs,  tm_a.*dis.red.*ARs,  tm_s.*ARs, ...  
                                    dis.red.*(1-dis.trv1).*ARs,  (1-dis.trv1).*ARs,  tm_a.*dis.red.*(1-dis.trv1).*ARs,  tm_s.*(1-dis.trv1).*ARs];
    J(1:ln,15*ln+1:16*ln)       = diag(dis.nu.*onesn);
    %Sh
    J(ln+1:2*ln,ln+1:2*ln)      = -FOI - diag(dis.hrv1.*onesn);
    J(ln+1:2*ln,5*ln+1:13*ln)   = -[dis.red.*ARsh,  ARsh,  tm_a.*dis.red.*ARsh,  tm_s.*ARsh, ...  
                                    dis.red.*(1-dis.trv1).*ARsh,  (1-dis.trv1).*ARsh,  tm_a.*dis.red.*(1-dis.trv1).*ARsh,  tm_s.*(1-dis.trv1).*ARsh];
    %Sv
    J(2*ln+1:3*ln,ln+1:2*ln)    = diag(dis.hrv1.*onesn);
    J(2*ln+1:3*ln,2*ln+1:3*ln)  = -(1-dis.scv1).*FOI - diag(dis.nuv1.*onesn);
    J(2*ln+1:3*ln,5*ln+1:13*ln) = -[dis.red.*ARsv,  ARsv,  tm_a.*dis.red.*ARsv,  tm_s.*ARsv, ...  
                                    dis.red.*(1-dis.trv1).*ARsv,  (1-dis.trv1).*ARsv,  tm_a.*dis.red.*(1-dis.trv1).*ARsv,  tm_s.*(1-dis.trv1).*ARsv];
    %E
    J(3*ln+1:4*ln,1:ln)         = FOI; 
    J(3*ln+1:4*ln,ln+1:2*ln)    = FOI;
    J(3*ln+1:4*ln,3*ln+1:4*ln)  = -diag((sig1+sig2+sig3+sig4).*onesn);
    J(3*ln+1:4*ln,5*ln+1:13*ln) = [dis.red.*(ARs+ARsh),  (ARs+ARsh),  tm_a.*dis.red.*(ARs+ARsh),  tm_s.*(ARs+ARsh), ...  
                                   dis.red.*(1-dis.trv1).*(ARs+ARsh),  (1-dis.trv1).*(ARs+ARsh),  tm_a.*dis.red.*(1-dis.trv1).*(ARs+ARsh),  tm_s.*(1-dis.trv1).*(ARs+ARsh)];
    %Ev
    J(4*ln+1:5*ln,2*ln+1:3*ln)  = (1-dis.scv1).*FOI;
    J(4*ln+1:5*ln,4*ln+1:5*ln)  = -diag((sig1+sig2+sig3+sig4).*onesn);
    J(4*ln+1:5*ln,5*ln+1:13*ln) = [dis.red.*ARsv,  ARsv,  tm_a.*dis.red.*ARsv,  tm_s.*ARsv, ...  
                                   dis.red.*(1-dis.trv1).*ARsv,  (1-dis.trv1).*ARsv,  tm_a.*dis.red.*(1-dis.trv1).*ARsv,  tm_s.*(1-dis.trv1).*ARsv];
    %Ina
    J(5*ln+1:6*ln,3*ln+1:4*ln)     = diag(sig1.*onesn);
    J(5*ln+1:6*ln,5*ln+1:6*ln)     = -diag(dis.g1.*onesn);
    %Ins
    J(6*ln+1:7*ln,3*ln+1:4*ln)     = diag(sig2.*onesn);
    J(6*ln+1:7*ln,6*ln+1:7*ln)     = -diag((g2+h).*onesn);
    %Isa
    J(7*ln+1:8*ln,3*ln+1:4*ln)     = diag(sig3.*onesn);
    J(7*ln+1:8*ln,7*ln+1:8*ln)     = -diag(dis.g1.*onesn);
    %Iss
    J(8*ln+1:9*ln,3*ln+1:4*ln)     = diag(sig4.*onesn);
    J(8*ln+1:9*ln,8*ln+1:9*ln)     = -diag((g2+h).*onesn);
    %Inav
    J(9*ln+1:10*ln,4*ln+1:5*ln)    = diag(sig1.*onesn);
    J(9*ln+1:10*ln,9*ln+1:10*ln)   = -diag(dis.g1.*onesn);
    %Insv
    J(10*ln+1:11*ln,4*ln+1:5*ln)   = diag(sig2.*onesn);
    J(10*ln+1:11*ln,10*ln+1:11*ln) = -diag((dis.g2_v1+dis.h_v1).*onesn);
    %Isav
    J(11*ln+1:12*ln,4*ln+1:5*ln)   = diag(sig3.*onesn);
    J(11*ln+1:12*ln,11*ln+1:12*ln) = -diag(dis.g1.*onesn);
    %Issv
    J(12*ln+1:13*ln,4*ln+1:5*ln)   = diag(sig4.*onesn);
    J(12*ln+1:13*ln,12*ln+1:13*ln) = -diag((dis.g2_v1+dis.h_v1).*onesn);
    %H
    J(13*ln+1:14*ln,6*ln+1:7*ln)   = diag(h.*onesn);
    J(13*ln+1:14*ln,8*ln+1:9*ln)   = diag(h.*onesn);
    J(13*ln+1:14*ln,13*ln+1:14*ln) = -diag((g3+mu).*onesn);
    %Hv
    J(14*ln+1:15*ln,10*ln+1:11*ln) = diag(dis.h_v1.*onesn);
    J(14*ln+1:15*ln,12*ln+1:13*ln) = diag(dis.h_v1.*onesn);
    J(14*ln+1:15*ln,14*ln+1:15*ln) = -diag((g3+mu).*onesn);
    %R
    J(15*ln+1:16*ln,5*ln+1:6*ln)   = diag(dis.g1.*onesn);
    J(15*ln+1:16*ln,6*ln+1:7*ln)   = diag(g2.*onesn);
    J(15*ln+1:16*ln,7*ln+1:8*ln)   = diag(dis.g1.*onesn);
    J(15*ln+1:16*ln,8*ln+1:9*ln)   = diag(g2.*onesn);
    J(15*ln+1:16*ln,13*ln+1:14*ln) = diag(g3.*onesn);
    J(15*ln+1:16*ln,15*ln+1:16*ln) = -diag(dis.nu.*onesn);
    %Rv
    J(16*ln+1:17*ln,9*ln+1:10*ln)  = diag(dis.g1.*onesn);
    J(16*ln+1:17*ln,10*ln+1:11*ln) = diag(dis.g2_v1.*onesn);
    J(16*ln+1:17*ln,11*ln+1:12*ln) = diag(dis.g1.*onesn);
    J(16*ln+1:17*ln,12*ln+1:13*ln) = diag(dis.g2_v1.*onesn);
    J(16*ln+1:17*ln,14*ln+1:15*ln) = diag(g3.*onesn);
    %D
    J(17*ln+1:18*ln,13*ln+1:14*ln) = diag(mu.*onesn);
    J(17*ln+1:18*ln,14*ln+1:15*ln) = diag(mu.*onesn);
    
    r = eigs(J,1,'largestreal');

end
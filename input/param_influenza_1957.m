function dis = param_influenza_1957

dis = struct;

%Probabilities
dis.ps  = 0.669;
dis.ihr = dis.ps*13.5*[0.0001 0.0001 0.0001 0.0001 ...
                       0.0001 0.0001 0.0001 0.0001 ...
                       0.0001 0.0025 0.0025 0.0025 ...
                       0.0025 0.0200 0.0200 0.0200 0.0200];    
dis.ifr = dis.ps*[0.0001 0.0001 0.0001 0.0001 ...
                  0.0001 0.0001 0.0001 0.0001 ...
                  0.0001 0.0025 0.0025 0.0025 ...
                  0.0025 0.0200 0.0200 0.0200 0.0200];   

%Durations
dis.Tlat  = 1.1;
dis.Tinc  = 1.4;
dis.Tay   = 2.5;
dis.Tsr   = 2.5;
dis.Tsh   = 2.5;
dis.Threc = 5.0;
dis.Thd   = 5.0;
dis.Ti    = 365;

%Transmission
dis.red  = 0.58;
%dis.R0  = 1.8000;
dis.beta = 0.062526278615573;

end
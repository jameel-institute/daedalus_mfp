function dis = param_covid_wildtype

dis = struct;

%Probabilities
dis.ps  = 0.595;
dis.ihr = [0.000016 0.000016 0.000408 0.000408 ...	
           0.010400 0.010400 0.034300 0.034300 ...	
           0.042500 0.042500 0.081600 0.081600 ...	
           0.118000 0.118000 0.166000 0.166000 0.184000];
dis.ifr = [0.000016 0.000016 0.000070 0.000070 ...
           0.000309 0.000309 0.000844 0.000844 ...
           0.001610 0.001610 0.005950 0.005950 ...
           0.019300 0.019300 0.042800 0.042800 0.078000];

%Durations
dis.Tlat  = 4.6;
dis.Tinc  = 5.8;
dis.Tay   = 2.1;
dis.Tsr   = 4.0;
dis.Tsh   = 4.0;
dis.Threc = 12.0;
dis.Thd   = 12.0;
dis.Ti    = Inf;

%Transmission
dis.red  = 0.58;
%dis.R0  = 2.8700;
dis.beta = 0.057977824645671;

end
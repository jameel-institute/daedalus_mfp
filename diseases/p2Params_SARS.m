function dis = p2Params_SARS

dis = struct;

%Probabilities
dis.ps  = 0.867;%https://pmc.ncbi.nlm.nih.gov/articles/PMC3371799/pdf/04-1165.pdf
dis.ihr = [0.0578 0.0578 0.0578 0.0578 ...	
           0.0816 0.0816 0.0816 0.0816 ...	
           0.3026 0.3026 0.3026 0.3026 ...	
           0.8670 0.8670 0.8670 0.8670 0.6018];
dis.ifr = dis.ps*[0.017 0.017 0.017 0.017 ...
                  0.024 0.024 0.024 0.024 ...
                  0.089 0.089 0.089 0.089 ...
                  0.255 0.255 0.255 0.255 0.177];%https://pmc.ncbi.nlm.nih.gov/articles/PMC7169690/pdf/TMI-14-21.pdf

%Durations
dis.Tlat  = 6.868132;%taken from epireview%4.6;
dis.Tinc  = 5.42;%taken from morgenstern%4.0;
dis.Tay   = 2.1;%assumed same as Covid
dis.Tsr   = 4.0;%assumed same as Covid
dis.Tsh   = 3.99 - (dis.Tlat - dis.Tinc);%taken from morgenstern%3.75;
dis.Threc = 23.5;%taken from epireview%26.5-3.75;
dis.Thd   = 35.9;%taken from epireview%23.7-3.75;
dis.Ti    = 365;

%Transmission
dis.red  = 0;%assumed%0.58;
%dis.R0  = 1.2;%https://pmc.ncbi.nlm.nih.gov/articles/PMC7134599/pdf/main.pdf%1.7500;%3.0000;
dis.beta = 0.036384519941551;

end
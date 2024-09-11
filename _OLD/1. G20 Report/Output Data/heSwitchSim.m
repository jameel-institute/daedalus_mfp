function [g,h]=heSwitchSim(inp1,inp2)

    load('UK36.mat','data');
    numInt=2;
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepSwineFlu(data,numInt,inp1,inp2);

    %%
    
    pr.sw=1;%switching on
    xoptim=[ones(length(data.G),1);data.xmin'];
    tvec=[-60  pr.Tm  NaN  365*5+1];
    
    
        data.wfhAv=[zeros(1,length(data.G));data.wfhAv];
        pr.betamod(2)=0.70;%min(1.30*pr.betamod(2),1);
        
            
        
        

        
        
        
    

        
    %%    
        
    [f,g,tvec]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);

    plotMultiOutd_pp(f,xoptim,tvec,data,pr);
    
    %%
    
    h(1)=f(end,5);%deaths
    h(1)=h(1)/(sum(data.Npop)/(10^5));%deaths (per 100k)
    
    fullx=[ones(length(data.G),1),repmat(reshape(xoptim,length(data.G),numInt),1,length(tvec))];
    fullx=fullx(:,1:length(tvec)-1);
    dgva=repmat((12/365)*data.obj,1,length(tvec)-1);
    
    h(2)=(tvec(end)-tvec(1))*sum((12/365)*data.obj)-(tvec(2:end)-tvec(1:end-1))*sum(fullx.*dgva,1)';%GDP loss ($, million)
    h(2)=h(2)/1000;%GDP loss ($, billion)%h(2)=round(100*h(2)/2/(12*sum(data.obj)),4);%GDP loss (% per annum)
    
    if f(end,3)>10, h=NaN(2,1); end
    
    h2=round(100*h(2)/sum(2*12*data.obj/1000),4);
    figHandles=findobj('Type','figure');
    figure(figHandles(2));
    title(['\textbf{Deaths:} ' num2str(h(1)) ' (per 100k) \big/ \textbf{GDP Loss:} ' num2str(h2) ' (\% per biennium)']);
    
end
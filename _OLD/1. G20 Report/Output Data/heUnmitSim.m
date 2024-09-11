function [g,h]=heUnmitSim(inp1,inp2)

    load('CN36.mat','data');
    numInt=2;
    
    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepSARS(data,numInt,inp1,inp2);

    %%
    
    pr.sw=1;%switching on
    xoptim=[ones(2*length(data.G),1)];
    tvec=[-60  pr.Tm  NaN  365*5+1];

    
        data.wfhAv=[zeros(3,length(data.G))];
        pr.betamod(2)=1;pr.betamod(3)=0.80;%min(1.30*pr.betamod(2),1);
        

        
        

        
        
        
    

        
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

        inds=sign(f(:,4)-pr.Hmax);
        inds(inds==0)=-1;
        inds=1+find(diff(inds));

        tpts=f(inds,1);tpts=[tpts;tvec(end)];ntp=length(tpts(2:2:end));
        intl=sum(tpts(2:2:end)-tpts(1:2:2*ntp));
        h(2)=h(2)+intl*dot(12/365*data.obj,1-min(data.unmit,1));
    
    h(2)=h(2)/1000;%GDP loss ($, billion)%h(2)=100*h(2)/2/(12*sum(data.obj));%GDP loss (% per annum)
    
    h(3)=f(end,10:13)*data.lyl;
   
    %if f(end,3)>10, h=NaN(2,1); end
    
    h2=round(100*h(2)/sum(2*12*data.obj/1000),4);
    figHandles=findobj('Type','figure');
    figure(figHandles(2));
    title(['\textbf{Deaths:} ' num2str(h(1)) ' (per 100k) \big/ \textbf{GDP Loss:} ' num2str(h2) ' (\% per biennium)']);
    
end
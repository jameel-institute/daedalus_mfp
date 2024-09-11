function [g,h]=heSingleSim(xoptim,inp1,inp2)

    load('UK36.mat','data');
    numInt=length(xoptim)/length(data.G);

    [pr,vx,NN,n,ntot,na,NNbar,NNrep,Dout,beta]=hePrepSwineFlu(data,numInt,inp1,inp2);
    
    %%
    
    pr.sw=0;%switching off
    
    tvec=[-60  linspace(pr.Tm,  365*5+1,  numInt+1)];
    
    if numInt==1 && all(xoptim==ones(length(data.G),1))
        data.wfhAv=[zeros(1,length(data.G));data.wfhAv(1,:)];
        pr.betamod(2)=min(1.30*pr.betamod(2),1);
        
    elseif numInt==2 && all(xoptim==[data.xmin';ones(length(data.G),1)])
        tvec(3)=vx.end;
        data.wfhAv=[zeros(1,length(data.G));flipud(data.wfhAv)];
        pr.betamod(3)=0.70;%min(1.30*pr.betamod(3),1);
        if not(vx.end<tvec(end))
            numInt=1;
            tvec(3)=[];
            xoptim=xoptim(1:length(data.G));
            data.wfhAv(end,:)=[];
            pr.betamod(end)=[];
            data.obj=NaN*data.obj;
        end
        
    else
        data.wfhAv=[zeros(1,length(data.G));repmat(data.wfhAv(2,:),numInt,1)];
        
    end
    
    %%
    
    [f,g,~]=heRunCovid19(pr,vx,n,ntot,na,NN,NNbar,NNrep,Dout,beta,xoptim,tvec,0,data);

    plotMultiOutd_pp(f,xoptim,tvec,data,pr);
    
    %%
    
    h(1)=f(end,5);%deaths
    h(1)=h(1)/(sum(data.Npop)/(10^5));%deaths (per 100k)
    
    fullx=[ones(length(data.G),1),reshape(xoptim,length(data.G),numInt)];
    
    dgva=repmat((12/365)*data.obj,1,length(tvec)-1);
    
    h(2)=(tvec(end)-tvec(1))*sum((12/365)*data.obj)-(tvec(2:end)-tvec(1:end-1))*sum(fullx.*dgva,1)';%GDP loss ($, million)
    if numInt==1 && all(xoptim==ones(length(data.G),1))
        inds=sign(f(:,4)-pr.Hmax);
        inds(inds==0)=-1;
        inds=1+find(diff(inds));

        tpts=f(inds,1);
        intl=sum(tpts(2:2:end)-tpts(1:2:end));
        h(2)=h(2)+intl*dot(12/365*data.obj,1-min(data.unmit,1));
    end
    h(2)=h(2)/1000;%GDP loss ($, billion)%h(2)=round(100*h(2)/2/(12*sum(data.obj)),4);%GDP loss (% per annum)
    
    if f(end,3)>10, h=NaN(2,1); end
    
    h2=round(100*h(2)/sum(2*12*data.obj/1000),4);
    figHandles=findobj('Type','figure');
    figure(figHandles(2));
    title(['\textbf{Deaths:} ' num2str(h(1)) ' (per 100k) \big/ \textbf{GDP Loss:} ' num2str(h2) ' (\% per biennium)']);

end
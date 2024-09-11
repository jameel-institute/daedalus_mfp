function h=p2Sim_RC(inp1,inp2,inp3)

    load(strcat(inp1,'.mat'),'data');
    numInt=2;
    
    if strcmp(inp2,'Swine Flu');
        [data,Dout,dis,p2]=p2_Flu2009(data,numInt,inp3);
    elseif strcmp(inp2,'Spanish Flu');
        [data,Dout,dis,p2]=p2_Flu1918(data,numInt,inp3);
    elseif strcmp(inp2,'Covid');
        [data,Dout,dis,p2]=p2_CovidWT(data,numInt,inp3);
    elseif strcmp(inp2,'SARS');
        [data,Dout,dis,p2]=p2_SARS(data,numInt,inp3);
    end
    
    %%
    
    p2.sw=1;%switching on
    xoptim=[ones(length(data.obj),1);data.x_econ(2,:)'];
    tvec=[-60  p2.Tm  NaN  365*5+1];
    p2.betamod(2)=0.70;%min(1.30*pr.betamod(2),1);

    %%    
        
    [tvec,f,g]=p2Run(data,Dout,tvec,dis,xoptim,p2);

    %p2Plot(data,tvec,f,xoptim,p2);
    
    %%
    
    h(1)=f(end,5);%deaths
    h(3)=h(1)/(sum(data.Npop)/(10^5));%deaths (per 100k)
    
    fullx=[ones(length(data.obj),1),repmat(reshape(xoptim,length(data.obj),numInt),1,length(tvec))];
    fullx=fullx(:,1:length(tvec)-1);
    dgva=repmat(data.obj,1,length(tvec)-1);
    
    h(2)=(tvec(end)-tvec(1))*sum(data.obj)-(tvec(2:end)-tvec(1:end-1))*sum(fullx.*dgva,1)';%GDP loss ($, million)
    %h(2)=round(h(2),2);
    h(4)=100*h(2)/sum(2*365*data.obj);
    %h(4)=round(h(4),2);
    %h(2)=h(2)/1000;%GDP loss ($, billion)%h(2)=round(100*h(2)/2/(12*sum(data.obj)),4);%GDP loss (% per annum)
    
    %if f(end,3)>10, h=NaN(2,1); end
    
    %h2=round(100*h(2)/sum(2*365*data.obj/1000),4);
    %figHandles=findobj('Type','figure');
    %figure(figHandles(2));
    %title(['\textbf{Deaths:} ' num2str(h(1)) ' (per 100k) \big/ \textbf{GDP Loss:} ' num2str(h(2)) ' (\% per biennium)']);
    %title(['\textbf{Deaths:} ' num2str(round(h(1))) ' \big/ \textbf{GDP Loss:} \$' num2str(round(h(2),2)) ' million']);
end
function [ROI1,ROI2,ROI3,EPVgain2] = p2ROI(inp1,inp2,inp3,inp4a,inp4b,p1,p2,p3,n)
    
    load(strcat(inp1,'.mat'),'data');
    costs    = data.pppf*sum(data.prepcost,1);%all P2 measures
    startu   = 1*costs(1:2:5);
    annual   = 1*costs(2:2:6);
    is       = str2double(inp4a(end));
    ie       = str2double(inp4b(end));
    startu   = sum(startu(is:ie-1));
    annual   = sum(annual(is:ie-1));
    r        = 0.03;             
    EPV_cost = startu+annual*sum(1./(1+r).^[0:n-1]);
    
    sec_la   = p2Sim(inp1,inp2,inp3,inp4a,1,1);
    sec_lb   = p2Sim(inp1,inp2,inp3,inp4b,1,1);
    gain     = (sec_la-sec_lb)/3;
    EPVgain1 = p1*gain*sum(1./(1+r).^[0:n-1]);
    EPVgain2 = p2*gain*sum(1./(1+r).^[0:n-1]);
    EPVgain3 = p3*gain*sum(1./(1+r).^[0:n-1]);

    ROI1     = 100*(EPVgain1-EPV_cost)/EPV_cost;
    ROI2     = 100*(EPVgain2-EPV_cost)/EPV_cost;
    ROI3     = 100*(EPVgain3-EPV_cost)/EPV_cost;
    
end
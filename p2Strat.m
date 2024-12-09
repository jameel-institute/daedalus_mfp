function data = p2Strat(data,inp3)
    
    lx = length(data.B);
    
    if strcmp(inp3,'Elimination');
        data.xoptim = [ones(1*lx,1);data.x_econ(:,2);data.x_elim(:,1);NaN(lx,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;NaN(1,lx);zeros(1,lx)];
        data.imand  = [2];
        data.inext  = [2,2,3,2,5];
    elseif strcmp(inp3,'Economic Closures');
        data.xoptim = [ones(2*lx,1);data.x_econ(:,2);data.x_econ(:,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
        data.imand  = [3];
        data.inext  = [2,3,3,4,5];
    elseif strcmp(inp3,'School Closures');
        data.xoptim = [ones(2*lx,1);data.x_schc(:,2);data.x_schc(:,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
        data.imand  = [3];
        data.inext  = [2,3,3,4,5];
    elseif strcmp(inp3,'No Closures');
        data.xoptim = [ones(5*lx,1)];
        data.hw     = [zeros(5,lx)];
        data.imand  = [NaN];
        data.inext  = [5];
    else
        error('Unknown Mitigation Strategy!');
    end

end
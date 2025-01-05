function data = dd_set_strategy(data,inp3);
    
    data.inp3 = inp3;
    lx        = length(data.B);
    
    if strcmp(inp3,'No Closures');
        data.xoptim = [ones(lx,1);NaN(3*lx,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);NaN(3,lx);zeros(1,lx)];
        data.imand  = [NaN];%configuration numbers with mandatory distancing (if any)
        data.inext  = [5,6];%index corresponds to event number in event function, values number of next configuration
    elseif strcmp(inp3,'School Closures');
        data.xoptim = [ones(2*lx,1);data.x_schc(:,2);data.x_schc(:,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
        data.imand  = [3];
        data.inext  = [2,3,3,4,5,6];
    elseif strcmp(inp3,'Economic Closures');
        data.xoptim = [ones(2*lx,1);data.x_econ(:,2);data.x_econ(:,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
        data.imand  = [3];
        data.inext  = [2,3,3,4,5,6];
    elseif strcmp(inp3,'Elimination');
        data.xoptim = [ones(lx,1);data.x_econ(:,2);data.x_elim(:,1);NaN(lx,1);ones(lx,1)];
        data.hw     = [zeros(1,lx);data.wfh;data.wfh;NaN(1,lx);zeros(1,lx)];
        data.imand  = [2];
        data.inext  = [2,3,2,5,6];
    else
        error('Unknown Mitigation Strategy!');
    end
    
end
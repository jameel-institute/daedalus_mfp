function data = dd_set_strategy(data,inp3);

ln        = length(data.NNs);
lx        = length(data.obj);
data.inp3 = inp3;
int       = 5;

if strcmp(inp3,'No Closures');
    data.xconf = [ones(1,lx);NaN(3,lx);ones(1,lx)];
    data.hw    = [zeros(1,lx);NaN(3,lx);zeros(1,lx)];
    data.imand = [NaN];%configuration numbers with mandatory distancing (if any)
    data.ev_fn = @dd_set_rules_noclosures;%event functions are defined in dd_set_strategy.m
    data.inext = [5,6];%index corresponds to event number in event function, value corresponds to number of next configuration
elseif strcmp(inp3,'School Closures');
    data.xconf = [ones(2,lx);data.x_schc(:,2)';data.x_schc(:,1)';ones(1,lx)];
    data.hw    = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
    data.imand = [3];
    data.ev_fn = @dd_set_rules_seclosures;
    data.inext = [2,3,3,4,5,6];
elseif strcmp(inp3,'Economic Closures');
    data.xconf = [ones(2,lx);data.x_econ(:,2)';data.x_econ(:,1)';ones(1,lx)];
    data.hw    = [zeros(1,lx);data.wfh;data.wfh;data.wfh;zeros(1,lx)];
    data.imand = [3];
    data.ev_fn = @dd_set_rules_seclosures;
    data.inext = [2,3,3,4,5,6];
elseif strcmp(inp3,'Elimination');
    data.xconf = [ones(1,lx);data.x_econ(:,2)';data.x_elim(:,1)';NaN(1,lx);ones(1,lx)];
    data.hw    = [zeros(1,lx);data.wfh;data.wfh;NaN(1,lx);zeros(1,lx)];
    data.imand = [2];
    data.ev_fn = @dd_set_rules_elimination;
    data.inext = [2,3,2,5,6];
else
    error('Unknown Mitigation Strategy!');
end

Dvec      = zeros(ln,ln,int);
for i = 1:int;
    Dvec(:,:,i) = dd_calc_contacts(data,data.xconf(i,:)',data.hw(i,:));
end
data.Dvec = Dvec;

end
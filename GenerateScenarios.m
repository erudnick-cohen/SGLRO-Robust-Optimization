function activeList = GenerateScenarios(...
    xBest,activeList, constraintList, uUB,uLB, ...
    uDomainConstrList,rScenFun)
%this scenario generation method uses an iterative optimization with the
%goal of making sure that max scenarios that violate additional constratins
%for the current solution lead to the generation of additional scenarios
%that only consider those additional constraints not considered in the
%original max. The purpose of this is to esnure that the max scenario
%generation does not contract the space of scenarios we may need to still
%sample to much. For this reason, this scenario generation adds the
%randomly sampled scenario last, as scenario reduction will cause similar
%issues with the max scenarios otherwise.

xOld = xBest;
%Change this parameter to sample multiple scenarios per iteration
remainingIterCount = 1;
%Generate a candidate infeasible scenarios
infeasList = [];
uq = randInterval(uLB,uUB);
while isempty(infeasList) && remainingIterCount > 0
    feasU = 1;
    infeasList = [];
    %sample a scenario
    uq = rScenFun();
    if feasU == 1        
        %go through the list of constraints and see if x is infeasible for
        %this scenario
        for m = 1:length(constraintList)
            fun = constraintList{m};
            fval = fun(xBest,uq);
            if(fval > 1e-6)
                m
                fun
                infeasList(end+1) = m;
            end
        end
        remainingIterCount = remainingIterCount - 1;
    end
end
%Check if the infeasible list is empty, just quit out of function if so
if(remainingIterCount <= 0 && isempty(infeasList))
    return
end





exitflag = 0;
xquery = xBest;

opts = optimoptions(@fmincon,'Algorithm','sqp');


%define objective function for maximizing L1 norm (sum) of violated
%constraints
otherobj = @(x,u) 0
for iii = infeasList
    cfun = constraintList{iii};
    otherobj = @(x,u) otherobj(x,u)+cfun(x,u);
end


[x, fval, exitflag] = fmincon(@(xv) -otherobj(xOld,xv),uq,[],[],[],[],...
    uLB(:)',uUB(:)',@(xv) conFunU(xv,xOld,infeasList,constraintList,uDomainConstrList),opts);

%if the optimizer failed(such as due to iterations), we just use the query
%which we know is a feasible solution
initinfeasList = infeasList;
if(exitflag >= 1)
    %add scenario
    activeList{end+1} = struct('x', xOld, 'u', x,'constr',initinfeasList);
    %make check for if there are new constraints violated, which we can use
    %to generated more scenarios.
    nInfList = [];
    for m = 1:length(constraintList)
        fun = constraintList{m};
        if(fun(xBest,x) > 1e-6)
            nInfList(end+1) = m;
        end
    end
    csSet = initinfeasList;
    %loop and make more new secnarios until we either fail or run out of
    %new scenarios to consider
    sDiff = setdiff(nInfList,csSet);
    while(~isempty(sDiff))
        disp('chain maximizing another scenario')
        infeasList = sDiff;
        otherobj = @(x,u) 0
        for iii = infeasList
            cfun = constraintList{iii};
            otherobj = @(x,u) otherobj(x,u)+cfun(x,u);
        end
        [x, fval, exitflag] = fmincon(@(xv) -otherobj(xOld,xv),x,[],[],[],[],...
            uLB(:)',uUB(:)',@(xv) conFunU(xv,xOld,infeasList,constraintList,uDomainConstrList),opts);
        if(exitflag >= 1)
            activeList{end+1} = struct('x', xOld, 'u', x,'constr',infeasList);
            nInfList = [];
            for m = 1:length(constraintList)
                fun = constraintList{m};
                if(fun(xBest,x) > 1e-6)
                    nInfList(end+1) = m;
                end
            end
            csSet = [csSet sDiff];
            sDiff = setdiff(nInfList,csSet);
        else
            sDiff = [];
        end
    end

else
    %if we don't fail, we want to still add the query, to get a better
    %robust feasible region by "double dipping" our infeasible scenario
    %generation
end

%note that we put the sampled scenario last, after all scenarios generated
%from it
activeList{end+1} = struct('x', xOld, 'u', uq,'constr',initinfeasList);



end
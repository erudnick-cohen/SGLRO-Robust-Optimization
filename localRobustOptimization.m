function [xD scenList] = localRobustOptimization(...
    objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList,scenList,uNom)
%Scenario Generation with Local Robust Optimization
%objfun: The objective function, a function of x
%constraintList: List (cell array) of constraints in format g(x,u) <= 0
%IC: Initial conditions
%xLB: Lower bounds on x
%xUB: Upper bounds on x
%uLB: Lower bounds on u
%uUB: Upper bounds on u
%uDomainConstrList: Constraints defining domain of uncertain parameters
%scenList: Set of initial scenarios provided
%rScenFun: Function which randomly samples a scenario
%uNom: Nominal scenario to used by local robust optimization step


conMult = 1;
xD = IC(:)';
fopts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
uopts = optimoptions(@fmincon,'Algorithm','sqp');
%fopts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
%uopts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed');
%fopts = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
%uopts = optimoptions(@fmincon,'Algorithm','active-set','Display','iter-detailed');
valid = 0;

%defines constraints imposed by scenarios
    function [c ceq] = conFunX(x)
        ceq = [];
        c = [];

        for j = 1:length(scenList)
            
            for i = scenList{j}.constr
                fun = constraintList{i};
                c(end+1) = conMult*fun(x,scenList{j}.u);
            end
        end
    end
%defines constraints on what a valid scenario is
    function [c ceq] = UCon(u)
        ceq = [];
        c = [];
        for i = 1:length(uDomainConstrList)
            fun = uDomainConstrList{i};
            c(end+1) = conMult*fun(u);
        end
    end

while valid == 0
    valid = 1;
    
    %attempt to maximize each constraint
    for k = 1:length(constraintList)
        cfun = constraintList{k};
        cVioObj = @(u) -cfun(xD,u);
        [newP,fval] = fmincon(cVioObj,uNom,[],[],[],[],uLB(:),uUB(:),@(u) UCon(u),uopts);
        
        %if the constraint is violated, add the scenario that was
        %generated
        if(-fval > 1e-6)
            valid = 0;
            scenList{end+1} = struct('x',xD,'u',newP,'constr',k);
        end
    end
    if(valid == 0)
        
        %solve scenario robust optimization problem
        xD = fmincon(@(x) objfun(x),IC(:)',[],[],[],[],...
            [xLB(:)'],[xUB(:)'],@(x) conFunX(x),fopts);
    end
end  
        


end
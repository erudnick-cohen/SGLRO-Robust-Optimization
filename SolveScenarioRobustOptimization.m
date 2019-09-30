function [xBest, cost,grad,hessout] = SolveScenarioRobustOptimization(...
        xBest, objfun, activeList, constraintList,xLB, xUB,IC)
%Scenario Robust Optimization
%xBest: The current best soluton
%objfun: The objective function, a function of x
%activeList: The input set of scenarios
%constraintList: List (cell array) of constraints in format g(x,u) <= 0
%xLB: Lower bounds on x
%xUB: Upper bounds on x
%IC: initial conditions

    conMult = 1e0;
    bestCost = inf;
    cost = bestCost;
    disp('Updating robust solution')
    %define constraints
    function [c, ceq, DC, DCeq] = conFunX(x)
        ceq = [];
        cconlen = 0;
        for qqq = 1:length(activeList)
            cconlen = cconlen+length(activeList{qqq}.constr);
        end
        
        c = zeros(cconlen,1);
        ind = 1;
        %scenario robustness feasibility   
        for j = 1:length(activeList)
            for i = activeList{j}.constr
                fun = constraintList{i};
                c(ind) = fun(x(1:length(xLB)),activeList{j}.u);
                ind = ind+1;
            end
        end
        c = c*conMult;
    end
    %perform two iterations of minimization for the current scenarios
    %one using xBest to account for new feasibility info, the other with a
    %random starting point to ensure globality of search
    opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');%,'MaxFunctionEvaluations',1e6,'MaxIterations',1e3);

    [x2, fval2, exitflag2,~,~,grad,hessout] = fmincon(@(x) objfun(x),[IC(:)'],[],[],[],[],...
        [xLB(:)'],[xUB(:)'],@(x) conFunX(x),opts);
    if(exitflag2 >0)
        if(bestCost < fval2)
        else
            xBest = x2(1:end);
            cost = fval2;
        end
    end

    
end
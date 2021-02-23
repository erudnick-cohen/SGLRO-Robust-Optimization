function [xBest, cost, remainingIterCount] = GetGoodFeasibleXOrDeath(...
        xBest, objfun, activeList, objScenList, constraintList,xLB, xUB,remainingIterCount,grads,IC,bestCost)
    %this function performs regular scenario based robust optimization
    disp('Updating robust solution')
    cost = bestCost;
    %define constraints
    function [c, ceq, DC, DCeq] = conFunX(x)
        ceq = [];
        c = zeros(length(activeList)*length(constraintList)+length(objScenList),1);
        ind = 1;
        %scenario robustness feasibility
        for i = 1:length(constraintList)
            for j = 1:length(activeList)
                fun = constraintList{i};
                c(ind) = fun(x(1:length(xLB)),activeList{j}.u);
                ind = ind+1;
            end
        end
        %Objective Scenario robustness for design
        for i = 1:length(objScenList)
            c(ind) = objfun(x(1:length(xLB)),objScenList{i}.u)-x(end);
            ind = ind+1;
        end
        %objective cost in each scenario
        if nargout > 2
            DC = zeros(length(activeList)*length(constraintList)+length(objScenList),length(xLB)+1);
            DCeq = [];
            ind = 1;
            %scenario robustness feasibility
            for i = 1:length(constraintList)
                for j = 1:length(activeList)
                    fun = grads{3}{i};
                    DC(ind,:) = [fun(x(1:length(xLB)),activeList{j}.u)' 0];
                    ind = ind+1;
                end
            end
            %Objective Scenario robustness for design
            for i = 1:length(objScenList)
                fun = grads{1};
                DC(ind,:) = [fun(x(1:length(xLB)),objScenList{i}.u)' -1];
                ind = ind+1;
            end
            DC = DC';
        end
    end
    %define objective
    function [f gradf] = objfunwgradX(x)
        f = x(end);
        xzero = zeros(size(xLB));
        if nargout > 1
            gradf = [xzero 1];
        end
    end
    %perform two iterations of minimization for the current scenarios
    %one using xBest to account for new feasibility info, the other with a
    %random starting point to ensure globality of search
    opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
    if ~isempty(grads)
        opts = optimoptions(opts,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
    end
    [x2, fval2, exitflag2] = fmincon(@(x) objfunwgradX(x),[IC(:)' 0],[],[],[],[],...
        [xLB(:)' -inf],[xUB(:)' inf],@(x) conFunX(x),opts);
    if(exitflag2 >0)
        if(bestCost < fval2)
        else
            xBest = x2(1:end-1);
            cost = fval2;
        end
    end
    %The design found here SHOULD be added to the scenario list, BUT ONLY
    %WHEN IT IS NON-EMPTY, as it is the worst case scenario of that set
    %which should be added
    if(length(objScenList)> 0)
        wcCost = -inf;
        wcU = [];
        for i = 1:length(objScenList)
            funcost = objfun(xBest,objScenList{i}.u);
            if(wcCost < funcost)
                wcCost = funcost;
                wcU = objScenList{i}.u;
            end
        end
        %objScenList{end+1} = struct('x', xBest, 'u', wcU,'cost',wcCost);
        global worstHist;
        global IterTracker;
        worstHist{end+1} = {IterTracker, xBest, wcU, wcU,wcCost};
    end        
end
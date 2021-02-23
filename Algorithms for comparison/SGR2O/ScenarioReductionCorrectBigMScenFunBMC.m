function [activeList, remainingIterCount] = ScenarioReductionCorrectBigMScenFunBMC(...
        activeList, constraintList, uUB,uLB,...
        uDomainConstrList, remainingIterCount,grads,rScenFun,compVal)
    %set the big M value
    M = 1e3;
    minVio = 1e-3;
    
    %build list of scenarios to check
    scenList = flip(1:length(activeList));
    %scenList = randperm(length(activeList));
    %note this isn't being used as the true query
    query = uLB;
    infeasList = {};
    nInfeas = 0;
    outOfScen = 0;
    while(length(infeasList) < 2 && remainingIterCount > 0)
        invalidScenario = 0;
        infeasList = {};
        nInfeas = 0;
        if(~isempty(scenList))
            query = activeList{scenList(end)}.u;
            scenList(end) = [];
        else
            outOfScen = 1;
            query = rScenFun();
            %make sure this query is feasible
            for k = 1:length(uDomainConstrList)
                fun = uDomainConstrList{k};
                if(fun(query)>compVal)
                    invalidScenario = 1;
                end
            end
        end
        if(invalidScenario == 0)
            for i = 1:length(activeList)
                xcheck = activeList{i}.x;
                violated = [];
                csum = 0;
                for j = 1:length(constraintList)
                    fun = constraintList{j};
                    funval = fun(xcheck,query);
                    if(funval > compVal)
                        violated(end+1) = j;
                        nInfeas = nInfeas+1;
                        csum = csum+funval;
                    end
                end
                if(~isempty(violated))
                    infeasList{end+1} = struct('scenario', i, 'constraints', violated, 'violationSum', csum, 'x', activeList{i}.x);
                end
            end
        end
        remainingIterCount = remainingIterCount - 1;
    end
    %quit if we ran out of time here
    if(remainingIterCount <= 0)
        return
    end
    %now that we have found a scenario that can replace our existing
    %scenarios, we refine it via an optimization problem
    %get list of scenarios to not use
    notList = [];
    for i = 1:length(infeasList)
        notList(end+1) = infeasList{i}.scenario;
    end
    objList = setdiff(1:length(activeList),notList);
    %first we define our new constraint function, which requires that all
    %infeasible conditions for this scenario remain the case and that the
    %scenario is feasible
        %define constraint function that includes scenarios and the
    %infeasibility constraints
    function [c, ceq, DC, DCeq] = conFun(x)
        lu = length(uLB);
        lo = length(objList);
        u = x(1:lu);
        z1 = x(lu+1:lu+lo);
        z2 = x(lu+lo+1:end);
        ceq = [];
        c = zeros(length(uDomainConstrList)+nInfeas+length(constraintList)*1*length(objList)+length(objList)*length(constraintList),1);
 
        for q = 1:length(uDomainConstrList)
            fun = uDomainConstrList{q};
            c(q) = fun(u);
        end
        %maintain infeasibility for scenarios being replaced
        ind = length(uDomainConstrList) + 1;
        for q = 1:length(infeasList)
            for v = infeasList{q}.constraints
                fun = constraintList{v};
                c(ind) = compVal - fun(activeList{infeasList{q}.scenario}.x,u);
                ind = ind+1;
            end
        end
        %define z2 for each scenario as being UPPER bounded by its
        %constraint value in its scenario
        %z2 has a maximum value of whatever our minimum constraint
        %violation being considered is
        indref = ind;
        for q = 1:length(objList)
            for v = 1:length(constraintList)
                fun = constraintList{v};
                %c(ind) = M*z(q) - fun(activeList{objList(q)}.x, u);
                c(ind) = -M*(1-z2(ind-indref+1)) + compVal - fun(activeList{objList(q)}.x, u);
                ind = ind+1;
            end
        end
        %define z1, which wants to be as large as possible, but has to be
        %less than M*z2
        indref = ind;
        for q = 1:length(objList)
            for v = 1:length(constraintList)
                c(ind) = z1(q)-z2(ind-indref+1);
                ind = ind+1;
            end
        end
        %         %define z2 for each scenario as being upper bounded by all
%         %contraints for that scenario
%         for q = 1:length(objList)
%             for v = 1:length(constraintList)
%                 fun = constraintList{v};
%                 %c(ind) = M*z(q) - fun(activeList{objList(q)}.x, u);
%                 c(ind) = M*z2(q) - fun(activeList{objList(q)}.x, u);
%                 ind = ind+1;
%             end
%         end
        %define the z logic constraint that only one z can be greater than
        %zero(currently not big M because I can't figure this one out for
        %some reason without having stupid numbers of contraints, so just
        %using the NLP cheating way
%         for q = 1:length(objList)
%             %c(ind) = (z2(q)-1/M^2*z2(q))+1/M*z1(q);
%             c(ind) = z1(q)*z2(q);
%             ind = ind+1;
%         end
%         %additionally, require that z1 must be LESS than z2
%         for q = 1:length(objList)
%             %c(ind) = (z2(q)-1/M^2*z2(q))+1/M*z1(q);
%             c(ind) = z1(q)-z2(q);
%             ind = ind+1;
%         end
        if nargout > 2
            DC = zeros(length(uDomainConstrList)+nInfeas+length(constraintList)*1*length(objList)+length(objList),...
                length(objList)+length(uLB)+length(objList)*length(constraintList));
            lu = length(uLB);
            DCeq = [];
            %xzeros = zeros(length(constraintList)*length(objList),1);
            for q = 1:length(uDomainConstrList)
                fun = grads{4}{q};
                DC(q,1:lu) = fun(u);
            end
            %maintain infeasibility for scenarios being replaced
            ind = length(uDomainConstrList) + 1;
            for q = 1:length(infeasList)
                for v = infeasList{q}.constraints
                    fun = grads{3}{v};
                    DC(ind,1:lu) = - fun(activeList{infeasList{q}.scenario}.x,u);
                    ind = ind+1;
                end
            end
            %define z2 for each scenario as being UPPER bounded by its
            %constraint value in its scenario
            %z2 has a maximum value of whatever our minimum constraint
            %violation being considered is
            indref = ind;
            for q = 1:length(objList)
                for v = 1:length(constraintList)
                     fun = grads{3}{v};
                    %c(ind) = M*z(q) - fun(activeList{objList(q)}.x, u);
                    DC(ind,1:lu) = -fun(activeList{objList(q)}.x, u);
                    DC(ind, lu+lo+ind-indref+1) = -M; 
                    ind = ind+1;
                end
            end
            %define z1, which wants to be as large as possible, but has to be
            %less than M*z2
            indref = ind;
            ind2 = 1;
            for q = 1:length(objList)
                for v = 1:length(constraintList)
                    fun = constraintList{v};
                    %c(ind) = M*z(q) - fun(activeList{objList(q)}.x, u);
                    DC(ind, lu+lo+ind2) = -1;
                    DC(ind,lu+q) = 1;
                    ind2 = ind2+1;
                end
                ind = ind+1;
            end
%             %define z for each scenario as being upper bounded by all
%             %contraints for that scenario
%             for q = 1:length(objList)
%                 for v = 1:length(constraintList)
%                     fun = grads{3}{v};
%                     %                     DC(ind,1:lu) = - fun(activeList{objList(q)}.x, u);
%                     %                     DC(ind,q) = M;
%                     DC(ind,1:lu) = -fun(activeList{objList(q)}.x, u);
%                     DC(ind,lu+length(objList)+q) = M;
%                     ind = ind+1;
%                 end
%             end
%             for q = 1:length(objList)
%                 DC(ind,:) = zeros(2*length(objList)+length(uLB),1);
%                 DC(ind,lu+q) = -z2(q);
%                 DC(ind,lu+length(objList)+q) = 1-z1(q);
%             end
%             for q = 1:length(objList)
%                 DC(ind,:) = zeros(2*length(objList)+length(uLB),1);
%                 DC(ind,lu+q) = 1;
%                 DC(ind,lu+length(objList)+q) = -1;
%             end
            DC = DC';
        end
        
    end

    %we next build our objective function, which is the sum of all
    %constraint violations in the currently non-violated scenarios. This objective is
    %defined to be maximization.
%     function out = objFun(u)
%         out = 0;
%         for q = 1:length(activeList)
%             for r = 1:length(constraintList)
%                 fun = constraintList{r};
%                 funval = fun(activeList{q}.x, u);
%                 out = out - funval*(funval>0);
%             end
%         end
%     end
    function [out, dout] = objFun(x)
        %out = -sum(x((length(uLB)+1):(length(uLB)+length(objList))));
        %out = -sum(x((length(uLB)+length(objList)+1):end));
        lu = length(uLB);
        lo = length(objList);
        u = x(1:lu);
        z1 = x(lu+1:lu+lo);
        z2 = x(lu+lo+1:end);
        out = -sum(z1);
        if nargout > 1
            dout = zeros(length(objList)+length(uLB)+length(objList)*length(constraintList),1);
            dout((lu+1):(lu+lo)) = -1;
            %dout((length(uLB)+1):end) = -1;
            %dout((length(uLB)+length(objList)+1):end) = -1;
        end
    end
        opts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','OptimalityTolerance',1e-3,'StepTolerance',1e-3);
        if ~isempty(grads)
            opts = optimoptions(opts,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
        end
        [sol, fval, exitflag] = fmincon(@(x) objFun(x),[query(:)' ones(length(objList),1)' zeros(length(objList)*length(constraintList),1)'],[],[],[],[],...
            [uLB(:)' zeros(length(objList),1)' zeros(length(objList)*length(constraintList),1)'],...
            [uUB(:)'  ones(length(objList),1)' ones(length(objList)*length(constraintList),1)'],@(x) conFun(x),opts);
        newQ = sol(1:length(uLB));
        %remove old scenarios from the active scenario list if refinement
        %suceeded
        if(exitflag < 1 && max(conFun(sol)) > compVal)
            %if we failed, just use the query
            newQ = query;
        else
            %if we didn't fail, then we add any scenarios with positive z
            %values to the deletion list
            z = sol((length(uLB)+1):end)
            for i = 1:length(objList)
                if(z(i)>0)
                    csum = 0;
                    violated = [];
                    for j = 1:length(constraintList)
                        fun = constraintList{j};
                        funval = fun(activeList{objList(i)}.x, newQ);
                        if(funval > minVio)
                            csum = csum+funval;
                            violated(end+1) = j;
                        end
                    end
                    if(csum >0)
                        infeasList{end+1} = struct('scenario', objList(i), 'constraints', violated, 'violationSum', csum, 'x', activeList{objList(i)}.x);
                    end
                end
            end
        end
    %also find the x with the lowest constraint violation
    minViolationScenario = 1;
    delList = [];
    for i = 1:length(infeasList)
        delList(end+1) = infeasList{i}.scenario;
        if(infeasList{i}.violationSum < infeasList{minViolationScenario}.violationSum)
            minViolationScenario = i;
        end
    end
    %if there is more than two scenarios, also find the x with the highest
%     %constraint violation
%     maxViolationScenario = 1;
%     for i = 1:length(infeasList)
%         if(infeasList{i}.violationSum > infeasList{maxViolationScenario}.violationSum)
%             maxViolationScenario = i;
%         end
%     end
    activeList(delList) = [];
    %add the new scenario
    activeList{end+1} = struct('x',infeasList{minViolationScenario}.x,'u',newQ);
    %if we have removed more than two scenarios, also keep the largest
    %violation x as an "extra" copy of the scenario in the active list
%     if(length(infeasList) > 2)
%         activeList{end+1} = struct('x',infeasList{maxViolationScenario}.x,'u',newQ);
%     end
    disp('reduced scenarios to:')
    newQ
end
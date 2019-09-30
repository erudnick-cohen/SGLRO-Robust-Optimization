function [xBest, activeList, gradOut,hessOut ] = scenarioGenerationRobustOptimization(...
    objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList, itercount,rScenFun)
%Scenario Generation Robust Optimization
%objfun: The objective function, a function of x
%constraintList: List (cell array) of constraints in format g(x,u) <= 0
%IC: initial conditions
%xLB: Lower bounds on x
%xUB: Upper bounds on x
%uLB: Lower bounds on u
%uUB: Upper bounds on u
%uDomainConstrList: Constraints defining domain of uncertain parameters
%itercount: Number of iterations to run for
%rScenFun: Function which randomly samples a scenario

remainingIterCount = itercount;


lenA = 0;
if(remainingIterCount > 0)
    activeList = {};
    xBest = IC;
    %Solve scenario robust optimization problem with no scenarios
    disp('Generating initial solution')
    [xBest, bCost,gradOut,hessOut] = SolveScenarioRobustOptimization(...
                xBest, objfun, {}, constraintList,xLB, xUB,xUB);
    while remainingIterCount > 0
        
        %find new infeasible scenario, handles both random sampling and
        %scenario generation
        disp('Searching for new infeasible scenario')
        activeList = GenerateScenarios(...
            xBest, activeList, constraintList, uUB,uLB,uDomainConstrList,rScenFun);
        
        if(lenA ~= length(activeList))
                %Solve scenario robust optimization problem with current
                %set of scenarios
                [xBest, bCost,gradOut,hessOut] = SolveScenarioRobustOptimization(...
                    xBest, objfun, activeList, constraintList,xLB, xUB,IC);
            disp('Current best robust cost seems to be')
            bCost
            %get lengths of each scenario list, to find out when a scenario
            %gets added
            lenA = length(activeList);
        end
        %Decrement iteration counter
        remainingIterCount = remainingIterCount - 1;
    end
else
    %if we could not find a feasible solution at all for any scenario
    xBest = nan;
    activeList = [];
end


end


function [xBest, activeList, objScenList ] = SGR20(...
    objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList,Nfeas,Nopt, itercount, feasCount, optCount,ScenDisc,grads,rScenFun,counters)
%objective function is the objective function as a function of x and u
%constraint list is a list of constraints as function handles of x and u
%LBs and UBs are bounds for each variable
%uDomainConstrList is the list of additional constraints defining the
%domain of possible scenarios
%N is the maximum number of scenarios to consider
%itercount is the number of iterations to run before stopping
remainingIterCount = itercount;

%define some global variables for performance tracking
global solHist;
solHist = {};
global infeasHist;
infeasHist = {};
global worstHist;
worstHist = {}

%transform objective function into an uncertain constraint by adding in one
%variable to the problem and one new contraint that defines that variable
%as the objective function, this needs to be done for computing 
disp('Generating initial solution')
% [xBest, uInit, remainingIterCount] = DetermineInitialSolutionWithObjective(...
%     objfun,constraintList,xLB,xUB, uUB,uLB, uDomainConstrList,...
%     remainingIterCount,grads);
lenS = 0;
lenA = 0;
%this iteration tracker isn't for the code, its for data analysis
global IterTracker;
IterTracker = 0;
if(remainingIterCount > 0)
    activeList = {};
    objScenList = {};
%     activeList{end+1} = struct('x',xBest,'u',uInit);
%     objScenList{end+1} = struct('x',xBest,'u',uInit,'cost',objfun(xBest,uInit));
    xBest = IC;
    objScenList{end+1} = struct('x', xBest, 'u', randInterval(uLB,uUB),'cost',10000);
    [xBest, bCost, remainingIterCount] = GetGoodFeasibleXOrDeathExtraIter(...
                xBest, objfun, {}, objScenList, constraintList,xLB, xUB,remainingIterCount,{},xUB,inf);
    while remainingIterCount > 0
        
        %find new infeasible scenario, do not consider objective function
        disp('Finding new infeasible solution')
        [xFeas, activeList, countLeft] = FindInfeasibleScenarioOnlyWithScenFunTripleDip(...
            xBest, activeList, constraintList,xLB,xUB, uUB,uLB,...
            uDomainConstrList, feasCount,grads,rScenFun);
        %remainingIterCount = remainingIterCount-feasCount+countLeft;

       %if there has been a change which would result in a new optimal
       %design, run normal scenario robust optimization to get a new
       %optimal solution
       %bCost = inf;

       %we choose xBest to be the best of the objective scenarios we have
       %post-pruning, as when the old xBest gets pruned, we need to stop
        %using it. If everything got pruned, then we just take xFeas as
        %xBest, as we know that it is feasible(or at least we think it is)
%         bCost = inf;
%         xBest = xFeas;
%         for i = 1:length(objScenList)
%             if(objScenList{i}.cost < bCost)
%                 xBest = objScenList{i}.x;
%             end
%         end
        
        if(lenA ~= length(activeList))
            [xBest, bCost, remainingIterCount] = GetGoodFeasibleXOrDeath(...
                xBest, objfun, activeList, objScenList, constraintList,xLB, xUB,remainingIterCount,grads,xBest,inf);
            [xBest, bCost, remainingIterCount] = GetGoodFeasibleXOrDeath(...
                xBest, objfun, activeList, objScenList, constraintList,xLB, xUB,remainingIterCount,grads,IC,bCost);
            disp('Current best robust cost seems to be')
            bCost
            %get lengths of each scenario list, to determine if we'll need to
            %do feasible search or not
            lenA = length(activeList);
        end
        
        counts = zeros(size(counters));
        for i = 1:length(counters)
            counts(i) = counters{i}.counter;
        end
        solHist{end+1} = {xBest,bCost, remainingIterCount,IterTracker,counts};
        %remainingIterCount = remainingIterCount-optCount+countLeft;
        %if we now have too many constraint scenarios, enter scenario reduction
        while(length(activeList) > Nfeas && remainingIterCount > 0)
            disp('Scenario Reduction: Feasibility')
            remainingIterCount
            [activeList, remainingIterCount] = ScenarioReductionCorrectBigMScenFunBMC(...
                activeList, constraintList, uUB,uLB,...
                uDomainConstrList, remainingIterCount*ScenDisc,grads,rScenFun,1e-6);
            remainingIterCount = remainingIterCount/ScenDisc;
        end
        
        if(lenA ~= length(activeList))
            [xBest, bCost, remainingIterCount] = GetGoodFeasibleXOrDeath(...
                xBest, objfun, activeList, objScenList, constraintList,xLB, xUB,remainingIterCount,grads,xBest,inf);
            [xBest, bCost, remainingIterCount] = GetGoodFeasibleXOrDeath(...
                xBest, objfun, activeList, objScenList, constraintList,xLB, xUB,remainingIterCount,grads,IC,bCost);
            disp('Current best robust cost seems to be')
            bCost
            %get lengths of each scenario list, to determine if we'll need to
            %do feasible search or not
            lenA = length(activeList);
        end
        IterTracker = itercount-remainingIterCount;
        remainingIterCount = remainingIterCount - 1;
    end
else
    %if we could not find a feasible solution at all for any scenario
    xBest = nan;
    activeList = [];
end


end


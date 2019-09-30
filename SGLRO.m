function [xRobust,uf] = SGLRO(objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList, itercount,rScenFun,uNom)
%Scenario Generation with Local Robust Optimization
%objfun: The objective function, a function of x
%constraintList: List (cell array) of constraints in format g(x,u) <= 0
%IC: Initial conditions
%xLB: Lower bounds on x
%xUB: Upper bounds on x
%uLB: Lower bounds on u
%uUB: Upper bounds on u
%uDomainConstrList: Constraints defining domain of uncertain parameters
%itercount: Number of iterations to run for
%rScenFun: Function which randomly samples a scenario
%uNom: Nominal scenario to used by local robust optimization step

%Sampling based robust optimization step
[xs1, uf1,~,~] = scenarioGenerationRobustOptimization(objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList, itercount,rScenFun);

%Local robust optimization step
[xRobust, uf] = localRobustOptimization(objfun,constraintList,xs1,xLB',xUB',uLB,uUB,uDomainConstrList,uf1,uNom);

end


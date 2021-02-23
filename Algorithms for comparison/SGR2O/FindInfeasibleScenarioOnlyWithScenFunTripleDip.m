function [xBest, activeList, remainingIterCount] = FindInfeasibleScenarioOnlyWithScenFunTripleDip(...
    xBest,activeList, constraintList,xLB, xUB, uUB,uLB, ...
    uDomainConstrList,remainingIterCount,grads,rScenFun)
xOld = xBest;
%Generate a candidate infeasible scenarios
infeasList = [];
uq = randInterval(uLB,uUB);
while isempty(infeasList) && remainingIterCount > 0
    feasU = 1;
    infeasList = [];
    uq = rScenFun();
    %check if uq is feasible(note that this should generally be the case
    %for this version, but this is here for more "heuristicish" approaches
    for m = 1:length(uDomainConstrList)
        fun = uDomainConstrList{m};
        if fun(uq) > 0
            feasU = 0;
        end
    end
    if feasU == 1        
        %go through the list of constraints and see if x is infeasible for
        %this
        for m = 1:length(constraintList)
            fun = constraintList{m};
            if(fun(xBest,uq) > 1e-6)
                infeasList(end+1) = m;
            end
        end
        remainingIterCount = remainingIterCount - 1;
    end
end
%also checking if the infeasible list is empty for sanity since we seemed
%to be passing this point when it way
if(remainingIterCount <= 0 && isempty(infeasList))
    return
end
%uq is the candidate infeasible scenario, we now solve the scenario
%minimization problem to find the scenario and design combo with the
%smalllest violation of these constraints for the OLD design(i.e. the
%scenario for which that design is closest to being feasible, but is
%still infeasible)

%define constraint function that includes scenarios and the
%infeasibility constraints
    function [c, ceq, DC, DCeq] = conFun(u)
        ceq = [];
        c = zeros(length(uDomainConstrList)...
            +length(infeasList),1);
        ind = 1;
        for i = 1:length(uDomainConstrList)
            fun = uDomainConstrList{i};
            c(ind) = fun(u);
            ind = ind+1;
        end
        %infeasibility
        for i = 1:length(infeasList)
            fun = constraintList{infeasList(i)};
            c(ind) = 1e-9 - fun(xOld,u);
            ind = ind+1;
        end
        if nargout > 2
            DC = zeros(length(uDomainConstrList)...
                +length(infeasList),length(uLB));
            DCeq = [];
            xzero = zeros(size(xLB));
            uzero = zeros(size(uLB));
            ind = 1;
            for i = 1:length(uDomainConstrList)
                fun = grads{7}(i);
                DC(ind,:) = fun(u);
                ind = ind+1;
            end
            %infeasibility
            for i = 1:length(infeasList)
                fun = grads{4}{infeasList(i)};
                DC(ind,:) = - fun(xOld,u);
                ind = ind+1;
            end
            DC = DC';
        end
    end
%define objective function of minimizing the infeasibility
    function [out dout] = conobjfun(u)
        ind = 1;
        out = 0;
        %infeasibility
        for i = 1:length(infeasList)
            fun = constraintList{infeasList(i)};
            out = out + fun(xOld,u);
            ind = ind+1;
        end
        if nargout >1
            dout = zeros(length(uLB),1);
            for i = 1:length(infeasList)
                fun = grads{4}{infeasList(i)};
                dout = dout + fun(xOld,u);
                ind = ind+1;
            end
        end
    end

%there is no optimization loop here because the initial point supplied is a
%feasible solution to the problem, if we fail to finish we just use it
exitflag = 0;
xquery = xBest;%randInterval(xLB,xUB);

opts = optimoptions(@fmincon,'Algorithm','sqp');
if ~isempty(grads)
    opts = optimoptions(opts,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
end
%with this implementation we "triple dip" we both min AND max the
%constraints.
%minimizing the contraints puts us on the feasible region, but it can move
%slowly because it depends on how close our x is to the robust feasible
%region
%using the sample preserves information
%maximizing speeds up the rate at which we push x to the boundary of the
%feasible region
[x, fval, exitflag] = fmincon(@(x) conobjfun(x),[uq],[],[],[],[],...
    uLB(:),uUB(:),@(x) conFun(x),opts);
[xM, fvalM, exitflagM] = fmincon(@(x) -conobjfun(x),[uq],[],[],[],[],...
    uLB(:),uUB(:),@(x) conFun(x),opts);
%if the optimizer failed(such as due to iterations), we just use the query
%which we know is a feasible solution
if(exitflag < 1)
    uInf = uq;
else
    %if we don't fail, we want to still add the query, to get a better
    %robust feasible region by "double dipping" our infeasible scenario
    %generation
    activeList{end+1} = struct('x', xOld, 'u', uq);
end
%if our max attempt succeeded, "triple dip" by using the constraint
%maximizing scenario as well
if(exitflagM > 0)
    activeList{end+1} = struct('x', xOld, 'u', xM);
end
uInf = x;
%note that we always put the min scenario last, as its the one that should
%be on the robust feasible region
activeList{end+1} = struct('x', xOld, 'u', uInf);
disp('new infeasible scenario')
activeList{end}
%     [x, fval, exitflag] = fmincon(@(x) optfun(x(1:length(xLB))),[x],[],[],[],[],...
%         [xLB uLB],[xUB uUB],@(x) conFun([x(1:length(xLB)) uInf]),opts);
%      [x, fval, exitflag] = fmincon(@(x) -optfun(x(1:length(xLB))),[x],[],[],[],[],...
%          [xLB uLB],[xUB uUB],@(x) conFun(x),opts);
xFeas = xBest
uq
x
xM
conFun(x)
infeasList
conFun(uq)
global infeasHist;
global IterTracker;
infeasHist{end+1} = {IterTracker, xOld,uq,conFun(uq),uInf,conFun(x)};


end
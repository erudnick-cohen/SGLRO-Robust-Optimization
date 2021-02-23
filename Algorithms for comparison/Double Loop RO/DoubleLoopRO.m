function [xD] = DoubleLoopRO(...
    objfun,constraintList,IC,xLB,xUB, uLB,uUB, uDomainConstrList)
conMult = 1;
xD = IC(:)';
pSet = ones(length(constraintList),length(uLB(:)')).*(uUB(:)'-uLB(:)')/2+uLB(:)';
pInit = pSet;
%fopts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
%uopts = optimoptions(@fmincon,'Algorithm','sqp','Display','iter-detailed');
fopts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e3);
uopts = optimoptions(@fmincon,'Algorithm','interior-point','Display','iter-detailed');
valid = 0;
    function [c ceq] = conFunX(x,pSet)
        ceq = [];
        c = zeros(length(constraintList),1);
        for i = 1:length(constraintList)
            fun = constraintList{i};
            c(i) = conMult*fun(x,pSet(i,:));
        end
    end
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

xD = fmincon(@(x) objfun(x),IC(:)',[],[],[],[],...
        [xLB(:)'],[xUB(:)'],@(x) conFunX(x,pSet),fopts);
    xD
    
    for j = 1:length(constraintList)
        cfun = constraintList{j};
        cVioObj = @(u) -cfun(xD,u);
        %[newP,fval] = fmincon(cVioObj,pSet(j,:),[],[],[],[],uLB(:),uUB(:),[],uopts);
        [newP,fval] = fmincon(cVioObj,pInit(j,:),[],[],[],[],uLB(:),uUB(:),@(u) UCon(u),uopts);
        newP
        %[newP,fval] = fmincon(cVioObj,pInit(j,:),[],[],[],[],uLB(:),uUB(:),[],uopts);
        %[newP,fval] = fmincon(cVioObj,randInterval(uLB,uUB),[],[],[],[],uLB(:),uUB(:),[],uopts);
        pSet(j,:) = newP;
        if(-fval > 1e-6)
            valid = 0;
        end
    end
end  
        


end
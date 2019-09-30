%define constraint function that includes scenarios and the
%infeasibility constraints
    function [c, ceq, DC, DCeq] = conFunU(u,xOld,CinfeasList,constraintList,uDomainConstrList)
        ceq = [];
        c = zeros(length(uDomainConstrList)...
            +length(CinfeasList),1);
        ind = 1;
        for i = 1:length(uDomainConstrList)
            fun = uDomainConstrList{i};
            c(ind) = fun(u(:));
            ind = ind+1;
        end
        %infeasibility
        for i = 1:length(CinfeasList)
            fun = constraintList{CinfeasList(i)};
            c(ind) = 1e-9 - fun(xOld(:),u(:));
            ind = ind+1;
        end
        if nargout > 2
            DC = zeros(length(uDomainConstrList)...
                +length(CinfeasList),length(uLB));
            DCeq = [];
            xzero = zeros(size(xLB));
            uzero = zeros(size(uLB));
            ind = 1;
            for i = 1:length(uDomainConstrList)
                fun = grads{7}(i);
                DC(ind,:) = fun(u(:));
                ind = ind+1;
            end
            %infeasibility
            for i = 1:length(CinfeasList)
                fun = grads{4}{CinfeasList(i)};
                DC(ind,:) = - fun(xOld(:),u(:));
                ind = ind+1;
            end
            DC = DC';
        end
    end
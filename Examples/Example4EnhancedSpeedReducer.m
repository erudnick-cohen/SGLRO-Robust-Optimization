clear x
xLB = [1 .1 8 0.3 0.3 1 1];
xUB = [15 1.5 28 12.3 12.3 8 8];

%uLB = [-.1 -0.01 -1 -.1 -.1 -.1 -.05];
uLB = [-.01 -0.01 -1 -.1 -.1 -.1 -0.05];
uUB = -uLB;
spaceTol = .3;

ufun = @(u) zeros(size(uLB));
scenFun = @() randInterval(uLB,uUB);
xu = @(x,u) (x+u(1:length(xLB)));
f0 = FunctionCounter();
f1 = FunctionCounter();
f2 = FunctionCounter();
f3 = FunctionCounter();
f4 = FunctionCounter();
f5 = FunctionCounter();
f6 = FunctionCounter();
f7 = FunctionCounter();
f8 = FunctionCounter();
f9 = FunctionCounter();
f10 = FunctionCounter();
f11 = FunctionCounter();
f12 = FunctionCounter();
f13 = FunctionCounter();
c1 = @(x) 1/(x(1)*x(2).^2.*x(3))-1/27;
c2 = @(x) 1/(x(1)*x(2).^2.*x(3).^2)-1/397.5;
c3 = @(x) x(4).^3/(x(2).*x(3).*x(6).^4)-1/1.93;
c4 = @(x) x(5).^3/(x(2).*x(3).*x(7).^4)-1/1.93;
c5 = @(x) x(2).*x(3)-40;
c6 = @(x) x(1)./x(2)-12;
c7 = @(x) 5-x(1)./x(2);
%c8 = @(x) (-x(4)+x(6))^2-.2;
%c9 = @(x) (-x(5)+x(7))^2-.2;
%c8 = @(x) 1.9-x(4)+1.5*x(6);
%c9 = @(x) 1.9-x(5)+1.1*x(7);
c8 = @(x) -x(4)+1*x(6)-spaceTol;
c9 = @(x) -x(5)+1*x(7)-spaceTol;
fo1 = @(x) .7864*x(1).*x(2).^2*(10/3*x(2).^2+14.933*x(3)-43.0934)-1.508*x(1)*(x(6).^2+x(7).^2)+7.477*(x(6).^3+x(7).^3)+.7864*(x(4).*x(6)^2+x(5).*x(7).^2);
fo2 = @(x) sqrt((745*x(4)./(x(2).*x(3))).^2+1.69e7)/(.1*x(6).^3);
fo3 = @(x) sqrt((745*x(5)./(x(2).*x(3))).^2+1.575e8)/(.1*x(7).^3);
c10 = @(x) fo2(x)-1800;
c11 = @(x) fo3(x)-1100;
c12 = @(x) fo1(x)-1400;
c13 = @(x,u) (x(4)-x(6)).^2+(x(5)-x(7)).^2-spaceTol^2;
radDist1 = @(x,u) (-u(8)*(x(4)-x(6)))/2;
radDist2 = @(x,u) (-u(9)*(x(5)-x(7)))/2;

%c12 = @(x,u) ((u(8)*(x(4)-x(6))/2)-u(9)*(x(5)-x(7))/2).^2-(spaceTol)^2;
cList = {};
cList{end+1} = @(x,u) f1.count(c1(xu(x,u)));
cList{end+1} = @(x,u) f2.count(c2(xu(x,u)));
cList{end+1} = @(x,u) f3.count(c3(xu(x,u)));
cList{end+1} = @(x,u) f4.count(c4(xu(x,u)));
cList{end+1} = @(x,u) f5.count(c5(xu(x,u)));
cList{end+1} = @(x,u) f6.count(c6(xu(x,u)));
cList{end+1} = @(x,u) f7.count(c7(xu(x,u)));
cList{end+1} = @(x,u) f8.count(c8(xu(x,u)));
cList{end+1} = @(x,u) f9.count(c9(xu(x,u)));
cList{end+1} = @(x,u) f10.count(c10(xu(x,u)));
cList{end+1} = @(x,u) f11.count(c11(xu(x,u)));
cList{end+1} = @(x,u) f12.count(c12(xu(x,u)));
cList{end+1} = @(x,u) f13.count(c13(xu(x,u),u));


uList = {};

%obj = @(xv,u) fo1(xu(xv,u));
%the new objective of this problem is to minimize the stresses present in
%the problem
obj = @(xv,u) (fo2(xu(xv,u))+fo3(xu(xv,u)))/1000; 
countList = {f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13};
ufun = @(u) zeros(size(u));
uNom = (uUB(:)'-uLB(:)')/2+uLB(:)';

xIC = [3.58, 0.71, 18, 8, 8, 3.5, 5.3];

dtimes = tic;
%Adjust the comments to change which algorithm you want to run

%Run SGLRO
[x, uf] = SGLRO(@(xv) f0.count(obj(xv,zeros(size(uLB)))), cList,xIC',xLB',xUB',uLB,uUB,uList,100,scenFun,uNom);

%Run Double Loop Approach (WARNING: WILL NOT TERMINATE)
%[x] = DoubleLoopRO(@(xv) f0.count(obj(xv,zeros(size(uLB)))),cList,xIC',xLB',xUB',uLB',uUB',uList)

%Run SGR20
%[x, uf,uo] = SGR20(@(xv,u) f0.count(obj(xv,ufun(u))), cList,xIC,xLB',xUB',uLB,uUB,uList,25,25,100,1,5,10,{},scenFun,countList)


runTime = toc(dtimes)

%%

%count function calls
splitObjectiveAndConstr
%%
%Verify feasibility and cost of final solution
objScenList = {};
objScenList{end+1} = struct('x', x, 'u', randInterval(uLB,uUB),'cost',10000);
uo = objScenList;
sList = {};
combList = {};
for i = 1:length(uLB)
    combList{end+1} = linspace(uLB(i),uUB(i),2);
end
allS = combvec(combList{:});
for i = 1:length(allS)
    sList{end+1} = struct('x', xLB, 'u', allS(:,i)');
end
uf1 = sList;
%%
worstCost = -inf;
for q = 1:length(uo)
    worstCost = max(worstCost,obj(x,ufun(uo{q}.u)));
end
worstCost

worstVio = -inf*ones(length(cList),1);
for q = 1:length(uf1)
    for j = 1:length(cList)
        fun = cList{j};
        worstVio(j) = max(worstVio(j),fun(x,uf1{q}.u));
    end
end
worstVio

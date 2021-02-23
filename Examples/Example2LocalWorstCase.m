clear x
xLB = [0];
xUB = [.8];

%uLB = [-.1 -0.01 -1 -.1 -.1 -.1 -.05];
uLB = [-1];
uUB = [1];

ufun = @(u) zeros(size(uLB));
scenFun = @() randInterval(uLB,uUB);
xu = @(x,u) (x+u);
f0 = FunctionCounter();
f1 = FunctionCounter();

c1 = @(x,u) (.75*u.*x.*(cos(7*pi/2*u.*(1-x))).^2-.1*x);

cList = {};
cList{end+1} = @(x,u) f1.count(c1(x,u));


uList = {};
obj = @(xv,u) -xv(1);%-xv(2).^2;
countList = {f0,f1};
ufun = @(u) zeros(size(u));

xIC = .1;
uNom = (uUB(:)'-uLB(:)')/2+uLB(:)';

dtimes = tic;
%Adjust the comments to change which algorithm you want to run

%Run SGLRO
[x, uf] = SGLRO(@(xv) f0.count(obj(xv,zeros(size(uLB)))), cList,xIC',xLB',xUB',uLB,uUB,uList,100,scenFun,uNom);

%Run Double Loop Approach 
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

clear x
nVar = 10; %Change this number to scale the problem sizes
nCon = nVar/2;% must be smaller than nVar
xLB = .1*ones(nVar,1);
xUB = .9*ones(nVar,1);

%uLB = [-.1 -0.01 -1 -.1 -.1 -.1 -.05];
uLB = -.09*ones(nVar,1);
uUB = -uLB;

ufun = @(u) zeros(size(uLB));
scenFun = @() randInterval(uLB,uUB);

g = @(xM) sum((xM(:)-.5).^2);

xu = @(x,u) (x+u);
f0 = FunctionCounter();
countList = {f0};
cList = {};
nm = nVar/nCon;

fBND = 1;
fj = @(x,j) sum(x((floor((j-1)*nm)+1):(floor((j)*nm))).^.9);
for i = 1:(nCon-1)
    f1 = FunctionCounter();
    nm = nVar/nCon;
    cfun  = @(x) fBND-fj(x,nCon)^2-fj(x,i)^2;
    cList{end+1} = @(x,u) f1.count(cfun(xu(x(:),u(:))));
    countList{end+1} = f1;
end
obj = @(x,u) sum(arrayfun(@(j) fj(x,j),1:nCon));

xIC = (xLB'+xUB')/2;
uNom = (uUB(:)'-uLB(:)')/2+uLB(:)';
scenFun = @() randInterval(uLB,uUB);

dtimes = tic;
%Adjust the comments to change which algorithm you want to run

%Run SGLRO
%[x, uf] = SGLRO(@(xv) f0.count(obj(xv,zeros(size(uLB)))), cList,xIC',xLB',xUB',uLB,uUB,uList,100,scenFun,uNom)

%Run Double Loop Approach
%[x] = DoubleLoopRO(@(xv) f0.count(obj(xv,zeros(size(uLB)))),cList,xIC',xLB',xUB',uLB',uUB',uList)

%Run SGR20
%[x, uf,uo] = SGR20(@(xv,u) f0.count(obj(xv,ufun(u))), cList,xIC,xLB',xUB',uLB,uUB,uList,25,25,100,1,5,10,{},scenFun,countList)


runTime = toc(dtimes)

%%

%count function calls
splitObjectiveAndConstr

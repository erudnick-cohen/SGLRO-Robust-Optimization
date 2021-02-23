clear x
%Define number of variables
%x = [h,l,t,b];
%define the uncertain parameters which are all deviations
%p = [h,l,t,b,F,L,taud,sigmad]
%define constants
F = 6000;
L = 14;
ct3 = .1047;
ct4 = .0481;
taud = 13600;
sigmad = 30000;
E = 30*10^6;
G = 12*10^6;
%define objective as summing all the x variables
obj = @(x,u) (1+ct3)*(x(1)+u(1)).^2.*(x(2)+u(2))+ct4*(x(3)+u(3)).*(x(4)+u(4))*(L+u(6)+x(2)+u(2));

cList = {};
%from that paper from 1976, which is referencing what I believe was the
%first edition of shigley
taup = @(x,u) ((F+u(5))/((sqrt(2)*(x(1)+u(1)).*(x(2)+u(2)))));
R = @(x,u) sqrt(1/4*(x(2)+u(2)).^2+((x(3)+u(3)+x(1)+u(1))/2).^2);
taupp = @(x,u) ((F+u(5))*((L+u(6)+(x(2)+u(2))/2)))*R(x,u)/...
    (2*(.707*(x(1)+u(1))*(x(2)+u(2))*(((x(2)+u(2)).^2)/12+(x(3)+u(3)+x(1)+u(1))/2).^2));
%c0 = @(x,u) obj(x(1:end-1),u)-x(end);
c0 = @(x,u) obj(x(1:end-1),u);
c1 = @(x,u) ((taup(x,u)^2+2*taup(x,u)*taupp(x,u)*(x(2)+u(2))/(2*R(x,u))+taupp(x,u)^2)^.5)/(taud+u(7)) - 1;
c2 = @(x,u) (6*(F+u(5))*(L+u(6))/((x(4)+u(4))*(x(3)+u(3))^2))/(sigmad+u(8)) - 1;
c3 = @(x,u) 4*(F+u(5))*(L+u(6))^3/(E*(x(3)+u(3))^3*(x(4)+u(4)))/.25 - 1;
I = @(x,u) 1/12*(x(3)+u(3))*(x(4)+u(4))^3;
alpha = @(x,u) 1/3*G*(x(3)+u(3))*(x(4)+u(4))^3;
c4 = @(x,u) (F+u(5))/(4.014*(sqrt(E*I(x,u)*alpha(x,u))/(L+u(6))^2)*(1-(x(3)+u(3))/(2*(L+u(6)))*sqrt(E*I(x,u)/alpha(x,u))))-1;
c5 = @(x,u)(x(1)+u(1))/(x(4)+u(4)) - 1;
c6 = @(x,u) .125/(x(1)+u(1)) - 1;
f0 = FunctionCounter();
f1 = FunctionCounter();
f2 = FunctionCounter();
f3 = FunctionCounter();
f4 = FunctionCounter();
f5 = FunctionCounter();
f6 = FunctionCounter();
%cList{end+1} = @(x,u) f0.count(c0(x,u));
cList{end+1} = @(x,u) f1.count(c1(x,u));
cList{end+1} = @(x,u) f2.count(c2(x,u));
cList{end+1} = @(x,u) f3.count(c3(x,u));
cList{end+1} = @(x,u) f4.count(c4(x,u));
cList{end+1} = @(x,u) f5.count(c5(x,u));
cList{end+1} = @(x,u) f6.count(c6(x,u));

uLB = [-.05 -.05 -.05 -.05 -500 -.2 -1000 -1000];
uUB = [.05 .05 .05 .05 500 .2 1000 1000];
ufun = @(u) zeros(size(u));


xLB = [.1 .1 .1 .1];
xUB = [2 10 10 2];
uList = {};
countList = {f0,f1,f2,f3,f4,f5,f6};


ufun = @(u) zeros(size(u));

xIC = xLB;
uNom = (uUB(:)'-uLB(:)')/2+uLB(:)';
scenFun = @() randInterval(uLB,uUB);

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

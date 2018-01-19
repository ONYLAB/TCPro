function optimizationdriver()

% RhoNT: maximum proliferation rate for activated helper T cells
RhoNT=0.1432; % day-1

% RhoAT: maximum proliferation rate for activated helper T cells
RhoAT=1.1; % day-1

%BetaAT: death rate of activated helper T cells
BetaAT=0.18; % day-1

A = [];
b = [];
Aeq = [];
beq = [];
lb = [0.1 0.5973 0.13];
ub = [0.2963 1.6 0.23];
nonlcon = [];
x0 = [RhoNT,RhoAT,BetaAT];

options = optimoptions('fmincon','Display','iter');

for i = 1:5
    x0 = fmincon(@(x) LamberthFig4B(x),x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    x0 = x0(:)'.*(1+(0.5-rand(1,3)));
end

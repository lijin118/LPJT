function [alpha,beta] = myQPSolver(Ys,Yt0,D,delta)

ns = size(Ys,1);
nt = size(Yt0,1);
labels = unique(Ys);
C = length(labels);

f = [];
A = [];
b = [];

Aeq_alpha = OneOfKEncoding(Ys)';
Aeq_beta = OneOfKEncoding(Yt0)';

Aeq_alpha = [Aeq_alpha;zeros(C-size(Aeq_alpha,1),ns)];
Aeq_beta = [Aeq_beta;zeros(C-size(Aeq_beta,1),nt)];
Aeq = [Aeq_alpha,zeros(C,nt);zeros(C,ns),Aeq_beta];

beq_alpha = delta*(sum(Aeq_alpha,2));
beq_beta = delta*(sum(Aeq_beta,2));

beq = [beq_alpha;beq_beta];

lb = zeros(ns+nt,1);
ub = ones(ns+nt,1);

options = optimset('Display','off', 'MaxIter', 1500, 'Algorithm', 'interior-point-convex','TolFun',1e-15);      
D = (D + D')*0.5;
[alpha_beta,~,exitflag] = quadprog(D,f,A,b,Aeq,beq,lb,ub,[],options);
alpha = alpha_beta(1:ns);
beta = alpha_beta(ns+1:end);
end
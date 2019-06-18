function [accVec] = LPJT_HDA(DataArr,options,silent)

Xs = DataArr.Xs;
Xt = DataArr.Xu;
Ys = DataArr.Ys;
Yt = DataArr.Yu;
Xl = DataArr.Xl;
Yl = DataArr.Yl;
ds = size(Xs,1);
dt = size(Xt,1);
ns = size(Xs,2);
nt = size(Xt,2);
nl = size(Xl,2);

T = options.T;
ReducedDim = options.ReducedDim;
mu = options.mu;
gamma = options.gamma;

delta = options.delta;
alpha = ones(ns,1);
beta = ones(nt,1);
accVec = zeros(1,T);

X = [Xl';Xt'];
lboptions= [];
lboptions.WeightMode = 'HeatKernel';
Wt = constructW(X,lboptions);
fl = zeros(nl,10);
labelIndex = (1:nl)';
fl(sub2ind(size(fl),labelIndex,Yl(labelIndex))) = 1;
[~, fu_CMN] = LabelPropagation(Wt, fl, 100);
[~,Yt0] = max(fu_CMN,[],2);

[Lws,Lbs] = myConGraph2(Ys,options,Xs');
Sbs = Xs*Lbs*Xs';
Sws = Xs*Lws*Xs';
Ht = eye(nt) - (1/nt)*ones(nt,nt);
Sht = Xt*Ht*Xt';

srcIndex = (1:ns)';
lblIndex = (1:nl)';
for i = 1:T
        [Lwt,Lbt] = myConGraph2(Yt0,options,Xt');

        Sbt = Xt*Lbt*Xt';
        Swt = Xt*Lwt*Xt';

        [g1,g2,Hss,Htt,Hst] = constructH(ns,nt,Ys,Yt0,delta);
        Mss = Xs*((alpha*alpha').*Hss)*Xs';
        Mtt = Xt*((beta*beta').*Htt)*Xt';
        Mst = -Xs*((alpha*beta').*Hst)*Xt'; 
        Mts = Mst';

        Smax = [gamma*Sbs,zeros(ds,dt);zeros(dt,ds),gamma*Sbt+mu*Sht];
        Smin = [Mss+gamma*Sws,Mst;Mts,Mtt+gamma*Swt+mu*eye(dt,dt)];

        [W,~] = eigs(Smax, Smin+eps*eye(ds+dt), ReducedDim, 'LM');   
        W = real(W);
        A = W(1:ds, :);
        B = W(ds+1:end, :);
        Zs = A'*Xs;
        Zt = B'*Xt;
        Zl = B'*Xl;

        X = [Zs';Zl';Zt'];

        Wt = constructW(X,lboptions);
        fl = zeros(ns+nl,10);
        fl(sub2ind(size(fl),srcIndex,Ys(srcIndex))) = 1;       
        fl(sub2ind(size(fl),ns+lblIndex,Yl(lblIndex))) = 1;        
        [~, fu_CMN] = LabelPropagation(Wt, fl, 100);
        [~,owner] = max(fu_CMN,[],2);
        Yt0 = owner;
        acc=length(find(Yt==Yt0))/length(Yt);
        accVec(i) = acc*100;
    
        if silent == 0
            fprintf('Accuracy of iteration %d is %.2f%% \n', i, accVec(i));      
        end

        G = [Xs'*A; Xt'*B] * [A'*Xs, B' * Xt];
        G1 = G.*g1;
        G2 = G.*g2;
        Kss = G1(1:ns,1:ns);
        Ktt  = G1(ns+1:end,ns+1:end);
        Kst = G2(1:ns, ns+1:end);
        Kts = Kst';
        D = [Kss,-Kst;-Kts,Ktt];
        delta = options.delta;  
        [alpha,beta] = myQPSolver(Ys,Yt0,D,delta);
end

end
function [accVec] = LPJT(DataArr,options)
Xs = DataArr.Xs;
Xt = DataArr.Xt;
Ys = DataArr.Ys;
Yt = DataArr.Yt;
m = size(Xs,1);
ns = size(Xs,2);
nt = size(Xt,2);

T = options.T;
ReducedDim = options.ReducedDim;
options.interK = floor(options.interK);
options.intraK = floor(options.intraK);
mu = options.mu;
lambda = options.lambda;
gamma = options.gamma;

delta = options.delta;
alpha = ones(ns,1);
beta = ones(nt,1);


[P,~] = pca1(Xs',options);
[Yt0,PCAAcc] = myClassifier(Xs,Xt,Ys,Yt,P);
fprintf('PCA accuracy : %f\n',PCAAcc);

accVec = zeros(T,1);
[Lws,Lbs] = myConGraph2(Ys,options,Xs');
Sbs = Xs*Lbs*Xs';
Sws = Xs*Lws*Xs';
Ht = eye(nt) - (1/nt)*ones(nt,nt);
Sht = Xt*Ht*Xt';

    for i = 1:T
        [Lwt,Lbt] = myConGraph2(Yt0,options,Xt');
        Sbt = Xt*Lbt*Xt';
        Swt = Xt*Lwt*Xt';       
        
        [g1,~,Hss,Htt,Hst] = constructH(ns,nt,Ys,Yt0,delta);

        Mss = Xs*((alpha*alpha').*Hss)*Xs';
        Mtt = Xt*((beta*beta').*Htt)*Xt';
        Mst = -Xs*((alpha*beta').*Hst)*Xt'; 
        Mts = Mst';      

        Smax = [gamma*Sbs,zeros(m,m);zeros(m,m),gamma*Sbt+mu*Sht];
        Smin = [Mss+gamma*Sws+lambda*eye(m),Mst-lambda*eye(m);Mts-lambda*eye(m),Mtt+gamma*Swt+(lambda+mu)*eye(m)];
        [W,~] = eigs(Smax, Smin+eps*eye(2*m), ReducedDim, 'LM');   
        A = W(1:m, :);
        B = W(m+1:end, :);
        Zs = A'*Xs;
        Zt = B'*Xt;               
        
        [accVec(i),Yt0] = mySVM(1,alpha,Zs',Ys,Zt',Yt,options.nu);
        fprintf('LPJT Accuracy of iteration %d: %.2f%% \n', i, accVec(i));                   
        
        G = [Xs'*A; Xt'*B] * [A'*Xs, B' * Xt];
        G1 = G.*g1;
        Kss = G1(1:ns,1:ns);
        Ktt  = G1(ns+1:end,ns+1:end);
        Kst = G1(1:ns, ns+1:end);
        Kts = Kst';
        D = [Kss,-Kst;-Kts,Ktt];
        [alpha,~] = myQPSolver(Ys,Yt0,D,delta);
        
    end

end

function [y,acc] = myClassifier(Xs,Xt,Ys,Yt,p)
        X = [Xs,Xt];
        Z = p'*X;
        Z = Z*diag(sparse(1./sqrt(sum(Z.^2))));
        Zs = Z(:,1:size(Xs,2));
        Zt = Z(:,size(Xs,2)+1:end);
        D=EuDist2(Zt',Zs');
         [~,idx]=sort(D,2);
         y=Ys(idx(:,1),1);
         acc=length(find(y==Yt))/length(Yt);
         clear X Z Zs Zt D Xs Xt Ys Yt p
end
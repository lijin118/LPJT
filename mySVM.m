function [acc,Yt0] = mySVM(scale,alpha,Zs,Ys,Zt,Yt,nu)

alpha = scale*alpha;

pivot_label = Ys;
pivot_feature = Zs;
pivot_weight = alpha;

svmoptions =strcat(['-q -s 1 -t 0 -n ',num2str(nu)]);
model = svmtrain(pivot_weight,pivot_label,sparse(pivot_feature),svmoptions);
[Yt0,acc,~] = svmpredict(Yt,sparse(Zt),model,'-q');
acc = acc(1);

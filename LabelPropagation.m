function [fu, fu_CMN] = LabelPropagation(W, fl, numIter)
% [fu, fu_CMN] = LabelPropagation(W, fl, numIter)
%
% Iterative version of LP (Zhu, 2002)
% Zhu, Xiaojin, and Zoubin Ghahramani. Learning from labeled and unlabeled data with label propagation. Technical Report CMU-CALD-02-107, Carnegie Mellon University, 2002.
%
% The parameters are equivalent to 'harmonic_function.m' by Zhu
% (And the additional parameter is numIter)
% numIter = forcing the algorithm is terminated on that iteration
%
% The code was implemented according to the pseudo code in Y. Bengio
% Y. Bengio et al., "Label propagation and quadratic criterion." Semi-supervised learning 10 (2006).
%
% Code ref: Fieping Nie, 'GeneralSSL.m'
% http://www.escience.cn/people/fpnie/papers.htm
%
% Byonghwa Oh, byonghwaoh@gmail.com, 2015

l = size(fl, 1); % the number of labeled points
n = size(W, 1); % total number of points
c = size(fl, 2); % the number of classes

D = sum(W)';
invD = spdiags(1./D, 0, n, n);

prevYu = [fl; zeros(n-l, c)];
Yu = zeros(n, c);

% Update iteratively
for i = 1:numIter    
    Yu = invD * W * prevYu;
    Yu(1:l, :) = fl;
    prevYu = Yu;
end

fu = Yu(l+1:n, :);

% compute the CMN solution (refer to the original harmonic function algorithm)
q = sum(fl) + 1;
fu_CMN = fu .* repmat(q./sum(fu), n-l, 1);

function [g1,g2,Hss,Htt,Hst] = constructH(ns,nt,Ys,Yt0,delta)

%%%%% Number of Categories %%%%%
lbl_idx = unique(Ys);
C = length(lbl_idx);

%%%% begin compute%%%%

% compute theta_m;
tm = ones(ns+nt,1);
tm(1:ns,1) = 1/(delta*ns);
tm(ns+1:end,1) = 1/(delta*nt);
theta_m = tm*tm';

theta_c1 = zeros(ns+nt,ns+nt);
% theta_c2 = zeros(ns+nt,ns+nt);
theta_c3 = zeros(ns+nt,ns+nt);
for i = 1:C
    Sc_index = Ys == lbl_idx(i);
    Tc_index = Yt0 == lbl_idx(i);
    nsc = sum(Sc_index);
    ntc = sum(Tc_index);
    if ntc == 0 || nsc ==0
        continue;
    end 
    
    ec = 1/(delta^2 * nsc * ntc);
    
    % theta_c1;
    theta = zeros(ns+nt,1);
    theta(Sc_index) = (1/(delta*nsc))*ones(nsc,1);
    theta(ns+find(Tc_index)) = (1/(delta*ntc))*ones(ntc,1);
    theta_c1 = theta_c1 + theta*theta';
    
    % theta_c3
    theta = zeros(ns+nt,1);
    theta(Sc_index) = ntc*ones(nsc,1);
    theta(ns+find(Tc_index)) = nsc*ones(ntc,1);
    theta_c3 = theta_c3 + diag(ec*(theta));
end

g1 = theta_m + theta_c1 + theta_c3;
g2 = theta_m + 2*theta_c1;
Hss = g1(1:ns,1:ns);
Htt = g1(ns+1:end,ns+1:end);
Hst = g2(1:ns, ns+1:end);
end
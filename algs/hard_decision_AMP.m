% This Matlab function generates the simulation results in the paper:
% 
% Jianan Bai and Erik G. Larsson, "Activity detection in distributed MIMO: 
% Distributed AMP via Likelihood Ratio Fusion,"
% in IEEE Wireless Communications Letters, 2022, doi: 10.1109/LWC.2022.3197053.
% 
% This is version 1.0 (Last edited: 2022-09-22)
% 
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% paper as described above.

function [decision,chEst_all,tau2_all] = hard_decision_AMP(Y_all,pilots_all,Rho,eps_all,maxIter,serviceMtx)

[L,KM] = size(Y_all);
[K,N_all] = size(serviceMtx);
M = KM/K;

lr_all = zeros(K,N_all);
chEst_all = zeros(N_all,KM);
tau2_all = zeros(K,1);
X_all = zeros(N_all,KM);
hard_decision = zeros(K,N_all);
B = zeros(K,N_all);
C = zeros(K,N_all);
for k = 1:K
    serviceList = find(serviceMtx(k,:));

    N = length(serviceList);

    pilots = pilots_all(:,serviceList);
    rho = Rho(k,serviceList);
    eps = eps_all(serviceList);
    
    X = zeros(N,M);
    Y = Y_all(:,(k-1)*M+1:k*M);
    Z = Y;
    
    tau2 = norm(Y,'fro')^2 / (L*M);
    tau2_record = [tau2];
    tau2_best = Inf;
    
    score_record = [];
    score_best = Inf;
    
    X_best = X;
    Z_best = Z;
    for iter = 1:maxIter
        Xi = X + pilots'*Z;
        Mtx = zeros(M,M);
        for n = 1:N
            xi = Xi(n,:).';
            omega = 1/tau2 - 1/(rho(n)+tau2);
            theta = 1/(1+(1-eps(n))/eps(n)*exp(M*log(rho(n)/tau2+1) - omega*norm(xi,2)^2));
            psi = rho(n)/(rho(n)+tau2);
            X(n,:) = theta*psi*xi;
            Mtx = Mtx + 1/L*theta*psi*(eye(M)+(1-theta)*omega*(xi*xi'));
        end
        Z = Y - pilots*X + Z*Mtx;
    
        tau2 = norm(Z,'fro')^2 / (L*M);
        tau2_record = [tau2_record,tau2];
        
        score = tau2;
        score_record = [score_record,score];
    
        if score < score_best
            score_best = score;
            tau2_best = tau2;
            X_best = X;
            Z_best = Z;
        end
    
        if score > 2*score_best
            break
        end
    end
    
    X = pilots'*Z_best + X_best;
    tau2 = tau2_best;
    tau2_all(k) = tau2;
    
    for n = 1:N
        x = X(n,:).';
        omega = 1/tau2 - 1/(rho(n)+tau2);
        lr = (rho(n)/tau2+1)^M*exp(-omega*norm(x,2)^2);
        lr_all(k,serviceList(n)) = lr;
        hard_decision(k,serviceList(n)) = (lr < eps(n)/(1-eps(n)));
        B(k,serviceList(n)) = tau2/rho(n)*log(1+rho(n)/tau2);
        C(k,serviceList(n)) = (1+tau2/rho(n))*log(1+rho(n)/tau2);
    end
    
    chEst_all(serviceList,(k-1)*M+1:k*M) = X_best;
    X_all(serviceList,(k-1)*M+1:k*M) = X;
end

Pmd = gammainc(M*ones(K,N_all),B*M,'lower')/gamma(M);
Pfa = gammainc(M*ones(K,N_all),C*M,'upper')/gamma(M);
Pmd = mean(Pmd,1).';
Pfa = mean(Pfa,1).';

thresh = (log((1-eps)./eps-M*log(Pmd./(1-Pfa)))) ./ log((1-Pmd).*(1-Pfa)./(Pmd.*Pfa));

decision = (sum(hard_decision,1).'>thresh);

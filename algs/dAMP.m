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

function [prob,chEst_all,tau2_all] = dAMP(Y_all,pilots_all,Rho,eps_all,maxIter,serviceMtx)

[L,KM] = size(Y_all);
[K,N_all] = size(serviceMtx);
M = KM/K;

prob_all = zeros(K,N_all);
chEst_all = zeros(N_all,KM);
tau2_all = zeros(K,1);
X_all = zeros(N_all,KM);
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
            gamma = 1/(1+(1-eps(n))/eps(n)*exp(M*log(rho(n)/tau2+1) - omega*norm(xi,2)^2));
            psi = rho(n)/(rho(n)+tau2);
            X(n,:) = gamma*psi*xi;
            Mtx = Mtx + 1/L*gamma*psi*(eye(M)+(1-gamma)*omega*(xi*xi'));
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
    
    chEst_all(serviceList,(k-1)*M+1:k*M) = X_best;
    X_all(serviceList,(k-1)*M+1:k*M) = X;
end

prob = zeros(N_all,1);
for n = 1:N_all
    APList = find(serviceMtx(:,n));
    AntInd = reshape((repmat((APList-1)*M+1,1,M) + (0:M-1)).',[],1);
    rho = Rho(APList,n);
    tau2t = tau2_all(APList);
    x = X_all(n,AntInd).';
    omega = 1./tau2t - 1./(rho+tau2t);
    omega_rep = reshape(repmat(omega,1,M).',[],1);
    prob(n) = 1/(1+exp(sum(M*log(rho./tau2t+1)) - norm(sqrt(omega_rep).*x,2)^2));
end
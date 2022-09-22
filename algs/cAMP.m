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

function [prob,chEst,tau2] = cAMP(Y,pilots,Rho,eps,maxIter,serviceMtx)

[L,KM] = size(Y);
[K,N] = size(serviceMtx);
M = KM/K;

X = zeros(N,KM);
Z = Y;

APList = cell(N,1);
for n = 1:N
    APList{n} = find(serviceMtx(:,n));
end

tau2 = zeros(K,1);
for k = 1:K
    tau2(k) = norm(Y(:,(k-1)*M+1:k*M),'fro')^2 / (L*M);
end
tau2_record = [tau2];
tau2_best = Inf(K,1);

score_record = [];
score_best = Inf(K,1);

for iter = 1:maxIter
    Xi = X + pilots'*Z;
    Mtx = zeros(KM,KM);
    for n = 1:N
        AntInd = reshape((repmat((APList{n}-1)*M+1,1,M) + (0:M-1)).',[],1);
        rho = Rho(APList{n},n);
        tau2t = tau2(APList{n});
        xi = Xi(n,AntInd).';
        omega = 1./tau2t - 1./(rho+tau2t);
        omega_rep = reshape(repmat(omega,1,M).',[],1);
        gamma = 1/(1 + (1-eps(n))/eps(n) * exp(...
            sum(M*log(rho./tau2t+1)) - norm(sqrt(omega_rep).*xi,2)^2));
        psi = rho ./ (rho+tau2t);
        psi_rep = reshape(repmat(psi,1,M).',[],1);
        X(n,AntInd) = gamma*psi_rep.*xi;
        Mtx(AntInd,AntInd) = Mtx(AntInd,AntInd) + 1/L*gamma*diag(psi_rep)*...
            (eye(length(AntInd)) + (1-gamma)*(xi*xi')*diag(omega_rep));
    end
    Z = Y - pilots*X + Z*Mtx;

    tau2 = zeros(K,1);
    for k = 1:K
        tau2(k) = norm(Z(:,(k-1)*M+1:k*M),'fro')^2 / (L*M);
    end
    tau2_record = [tau2_record,tau2];
    
    score = norm(Z,'fro')^2 / (L*KM);
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

prob = zeros(N,1);
for n = 1:N
    AntInd = reshape((repmat((APList{n}-1)*M+1,1,M) + (0:M-1)).',[],1);
    rho = Rho(APList{n},n);
    tau2t = tau2(APList{n});
    x = X(n,AntInd).';
    omega = 1./tau2t - 1./(rho+tau2t);
    omega_rep = reshape(repmat(omega,1,M).',[],1);
    prob(n) = 1/(1+exp(sum(M*log(rho./tau2t+1)) - norm(sqrt(omega_rep).*x,2)^2));
end

chEst = X_best;
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

function serviceMtx = GenerateServiceMtx_AMP(Rho,N_UE,N_AP,opt)

[K,N] = size(Rho);

serviceMtx = zeros(K,N);

if strcmp(opt,'AP')
    for k = 1:K
        [~,ind] = sort(Rho(k,:),'descend');
        serviceMtx(k,ind(1:N_UE)) = 1;
    end
elseif strcmp(opt,'UE')
    for n = 1:N
        [~,ind] = sort(Rho(:,n),'descend');
        serviceMtx(ind(1:N_AP),n) = 1;
    end
else
    error('Please specify the AP association option for AMP.')
end

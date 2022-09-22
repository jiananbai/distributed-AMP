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

function [serviceMtx,UEwoAP] = GenerateServiceMtx(nLSFC,txPower,SNR_th)

SNR = 10*log10(nLSFC.*txPower);

serviceMtx = (SNR > SNR_th);

UEwoAP = find(sum(serviceMtx,1)==0);
for i = UEwoAP
    [~,i_AP] = max(nLSFC(:,i));
    serviceMtx(i_AP,i) = 1;
end
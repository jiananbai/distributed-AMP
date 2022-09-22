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

function txPower = CalculateTxPower(nLSFC,maxPower,serviceMtx,UEwoAP,opt)

N = size(nLSFC,2);
UEwAP = 1:N;
UEwAP(UEwoAP) = [];

txPower = maxPower * ones(N,1);

s = zeros(N,1);
for n = UEwAP
    temp = nLSFC(serviceMtx(:,n)>0,n);
    if strcmp(opt,'FullPower')
        s(n) = 1;
    elseif strcmp(opt,'MasterAP')
        s(n) = max(temp);
    elseif strcmp(opt,'WorstAP')
        s(n) = min(temp);
    elseif strcmp(opt,'AvgAP')
        s(n) = mean(temp);
    elseif strcmp(opt,'SumAP')
        s(n) = sum(temp);
    else
        error('Please specify a power allocation option.')
    end

    txPower(UEwAP) = min(s(UEwAP))./s(UEwAP)*maxPower;
end
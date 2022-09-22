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

function dist = CalculateDistance(locAP,locUE,D)

K = size(locAP,1);
N = size(locUE,1);

clocAP = locAP(:,1) + 1j*locAP(:,2);
clocUE = locUE(:,1) + 1j*locUE(:,2);

trans_wrap = [-1;0;1] * D;
trans_wrap = reshape(trans_wrap + 1j*trans_wrap.',[],1);
clocAP_wrap = clocAP + trans_wrap.';

dist = zeros(K,N,9);
for i = 1:9
    dist(:,:,i) = abs(clocAP_wrap(:,i) - clocUE.');
end
dist = min(dist,[],3);
dist = sqrt(dist.^2 + 1e-4);  % the height of AP is 10 meters
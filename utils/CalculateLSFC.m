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

function LSFC = CalculateLSFC(dist,std_sf)

[K,N] = size(dist);

LSFC_dB = -140.6 - 36.7*log10(dist) + std_sf*randn(K,N);

LSFC = 10.^(LSFC_dB/10);
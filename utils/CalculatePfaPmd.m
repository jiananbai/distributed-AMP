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

function PfaPmd = CalculatePfaPmd(activity,stat,N_th)

thresh = linspace(0,max(stat,[],'all'),N_th);

N_active = sum(activity,'all');
[N_trials,N] = size(activity);
N_total = N_trials*N;

N_FA = zeros(N_th,1);
N_MD = zeros(N_th,1);
for i_th = 1:N_th
    decision = (stat > thresh(i_th));
    N_FA(i_th) = sum((decision-activity)>0,'all');
    N_MD(i_th) = sum((decision-activity)<0,'all');
end
P_FA = N_FA / (N_total-N_active);
P_MD = N_MD / N_active;

PfaPmd = [P_FA,P_MD];
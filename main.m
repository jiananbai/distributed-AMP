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

function main(L,powerOption)
% powerOption: 'FullPower', 'MasterAP', 'WorstAP', 'AvgAP', 'SumAP'

if nargin < 2
    L = 40;
    powerOption = 'FullPower';
end

addpath utils
addpath algs
addpath algs/ML

%% Parameter Settings
K = 20;  % number of APs
M = 3;  % number of antennas per AP
N = 400;  % number of devices
eps_0 = 0.1;  % average access probability
D = 2;  % area size in km
SNR_th = 6;  % SNR threshold for connecting to APs
maxPower_dBm = 23;  % maximum transmit power in dBm
maxPower = 10^(maxPower_dBm/10) * 1e-3;  % maximum transmit power in W
B = 1e6;  % transmit bandwidth in Hz
noise_PSD_dBm = -169;  % noise PSD in dBm/Hz
noise_power = 10^(noise_PSD_dBm/10) * B * 1e-3;  % noise power in W
std_sf = 4;  % shadow fading standard deviation in dB
maxIter = 10;  % maximum number of iterations for AD algorithms
N_th = 10000;  % number of threshold values for decision
N_AP_AMP = 10;  % number of APs to detect a UE in AMP
N_UE_AMP = 200;  % number of UEs that an AP will detect in AMP
 
%% Monte Carlo Simulation
N_trials = 5e5;

[gamma_record,probD1_record,probD2_record,probC1_record,probC2_record] =...
    deal(zeros(N_trials,N),zeros(N_trials,N),zeros(N_trials,N),zeros(N_trials,N),zeros(N_trials,N));

activity_record = zeros(N_trials,N);

[runtime_ML,runtime_D1,runtime_D2,runtime_C1,runtime_C2] = deal(0,0,0,0,0);

% parpool('local',parcluster('local').NumWorkers);
parfor trial = 1:N_trials
    eps = eps_0 * ones(N,1);

    pilots = GeneratePilots(L,N); 

    locAP = unifrnd(-D/2,D/2,K,2);
    locUE = unifrnd(-D/2,D/2,N,2);

    dist = CalculateDistance(locAP,locUE,D);

    LSFC = CalculateLSFC(dist,std_sf);
    nLSFC = LSFC/noise_power;

    [serviceMtx,UEwoAP] = GenerateServiceMtx(nLSFC,maxPower,SNR_th);
    txPower = CalculateTxPower(nLSFC,maxPower,serviceMtx,UEwoAP,powerOption);

    Rho = nLSFC * diag(txPower) * L;

    activeUE = binornd(ones(N,1),eps);
    supp = find(activeUE);
    activity_record(trial,:) = activeUE;

    H_all = zeros(N,K*M);
    Y_all = zeros(L,K*M);
    for k = 1:K
        H = sqrt(1/2) * (randn(N,M) + 1j*randn(N,M));
        H_all(:,(k-1)*M+1:k*M) = H;
        
        rho = Rho(k,:).';

        X = zeros(N,M);
        X(supp,:) = sqrt(diag(rho(supp))) * H(supp,:);

        W = sqrt(1/2) * (randn(L,M) + 1j*randn(L,M));
        
        Y = pilots * X + W;
        Y_all(:,(k-1)*M+1:k*M) = Y;

    end

    serviceMtx_AMP = GenerateServiceMtx_AMP(Rho,N_UE_AMP,N_AP_AMP,'UE');

%     [decision,~,~] = hard_decision_AMP(Y_all,pilots,Rho,eps,maxIter,ones(K,N));
    
    tic;
    gamma = decode_activity_pattern_Mpoly(1,Rho.',reshape(Y_all,L,M,K),pilots,3);
    runtime_ML = runtime_ML + toc/N_trials;

    tic; 
    [probC1,~,~] = cAMP(Y_all,pilots,Rho,eps,maxIter,ones(K,N));
    runtime_C1 = runtime_C1 + toc/N_trials;

    tic;
    [probD1,~,~] = dAMP(Y_all,pilots,Rho,eps,maxIter,ones(K,N));
    runtime_D1 = runtime_D1 + toc/N_trials;
    
    tic;
    [probC2,~,~] = cAMP(Y_all,pilots,Rho,eps,maxIter,serviceMtx_AMP);
    runtime_C2 = runtime_C2 + toc/N_trials;
    
    tic;
    [probD2,~,~] = dAMP(Y_all,pilots,Rho,eps,maxIter,serviceMtx_AMP);
    runtime_D2 = runtime_D2 + toc/N_trials;
    
    gamma_record(trial,:) = gamma;
    probD1_record(trial,:) = probD1;
    probD2_record(trial,:) = probD2;
    probC1_record(trial,:) = probC1;
    probC2_record(trial,:) = probC2;
    
end

PfaPmd_ML = CalculatePfaPmd(activity_record,gamma_record,N_th);
PfaPmd_D1 = CalculatePfaPmd(activity_record,probD1_record,N_th);
PfaPmd_D2 = CalculatePfaPmd(activity_record,probD2_record,N_th);
PfaPmd_C1 = CalculatePfaPmd(activity_record,probC1_record,N_th);
PfaPmd_C2 = CalculatePfaPmd(activity_record,probC2_record,N_th);

% result_path = append('results/result_L',num2str(L),'_',powerOption,'.mat');
% save(result_path, ...
%     'PfaPmd_ML', 'PfaPmd_D1','PfaPmd_D2','PfaPmd_C1','PfaPmd_C2', ...
%     'runtime_ML','runtime_D1','runtime_D2','runtime_C1','runtime_C2')
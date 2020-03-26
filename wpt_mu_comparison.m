clear; close all; clc; initialize; config_mu_comparison;
%% Waveform design by Max-Min-Rand, CHE Max-Min-Rand, and CHE Max-Min-RR algorithms
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
% [~, ~, ~, minVoltageRr] = waveform_max_min_rr(beta2, beta4, powerBudget, channel, tolerance, weight);
[~, ~, ~, minVoltageCheRr] = waveform_max_min_che_rr(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
minVoltageCheRand = zeros(nRealizations, 1);
for iRealization = 1 : nRealizations
    [~, ~, ~, minVoltageCheRand(iRealization)] = waveform_max_min_che_rand(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
end

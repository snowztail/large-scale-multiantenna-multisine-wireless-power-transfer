clear; close all; clc; initialize; config_max_min_rr;
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
nCandidates = 1e3;
[waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_rand(beta2, beta4, powerBudget, channel, tolerance, weight, nCandidates);

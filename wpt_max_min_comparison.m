clear; close all; clc; initialize; config_max_min_comparison;
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
[~, ~, ~, minVoltage] = waveform_max_min_che_rr(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);

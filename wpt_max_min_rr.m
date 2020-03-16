clear; close all; clc; initialize; config_max_min_rr;
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
[waveform, sumVoltage, userVoltage] = waveform_max_min_rr(beta2, beta4, powerBudget, channel, tolerance, weight);

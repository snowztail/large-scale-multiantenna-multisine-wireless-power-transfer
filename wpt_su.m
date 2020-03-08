clear; close all; clc;
%% Transceiver
powerBudget = db2pow(36 - 30);
distance = 10;
nSubbands = 8;
nTxs = 4;
nRxs = 1;
%% Channel
centerFrequency = 2.4e9;
bandwidth = 1e7;
carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
fadingType = "flat";
%% Harvester
% assumptions: antenna impedance = 50 ohms, ideality factor = 1, thermal voltage = 25.85 mV
beta2 = 9.6712e2;
beta4 = 6.0304e6;
tolerance = 1e-3;

% \boldsymbol{h}_{q,n}
subchannel = channel_tgn_e(distance, nSubbands, nRxs, nTxs, carrierFrequency, fadingType);
% s_n
waveform = waveform_su(beta2, beta4, powerBudget, subchannel, tolerance);
% v_{\text{out},q}
voltage = harvester(beta2, beta4, waveform, subchannel);

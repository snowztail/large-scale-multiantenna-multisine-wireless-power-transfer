clear; close all; clc; initialize; config_wsum;
%% Waveform design by WSum algorithm
channel = zeros(nTxs, nSubbands, nUsers);
for iUser = 1 : nUsers
    % \boldsymbol{h}_{q,n}
    channel(:, :, iUser) = channel_tgn_e(pathloss, nSubbands, nTxs, carrierFrequency, fadingType);
end
% \boldsymbol{s_n}
% waveform = waveform_wsum(beta2, beta4, powerBudget, channel, tolerance, weight);
% waveform1 = waveform_su(beta2, beta4, powerBudget, channel, tolerance);
waveform2 = waveform_wsums(beta2, beta4, powerBudget, channel, tolerance, weight);
[waveform3, asymWaveform3] = waveform_che_wsum(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);

% v_{\text{out},q}
% voltage = harvester(beta2, beta4, waveform, channel);
% voltage1 = harvester(beta2, beta4, waveform1, channel);
voltage2 = harvester_wsum(beta2, beta4, waveform2, channel, weight);
voltage3 = harvester_wsum(beta2, beta4, asymWaveform3, channel, weight);


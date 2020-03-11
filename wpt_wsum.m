clear; close all; clc; initialize; config_wsum;
%% Waveform design by WSum algorithm
channel = zeros(nTxs, nSubbands, nUsers);
for iUser = 1 : nUsers
    % \boldsymbol{h}_{q,n}
    channel(:, :, iUser) = channel_tgn_e(distance, nSubbands, nTxs, carrierFrequency, fadingType);
end
% \boldsymbol{s_n}
[waveform] = waveform_wsum(beta2, beta4, powerBudget, channel, tolerance, weight);
waveform1 = waveform_su(beta2, beta4, powerBudget, channel, tolerance);
% v_{\text{out},q}
voltage = harvester(beta2, beta4, waveform, channel);
voltage1 = harvester(beta2, beta4, waveform1, channel);

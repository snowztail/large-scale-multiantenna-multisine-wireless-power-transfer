clear; close all; clc; initialize; config_wsum;
%% Waveform design by SU WPT and WSum algorithms
voltageSu = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
voltageWsum = zeros(length(Variable.nTxs), length(Variable.nSubbands), nRealizations);
for iTx = 1 : length(Variable.nTxs)
    nTxs = Variable.nTxs(iTx);
    powerBudget = eirp / nTxs;
    for iSubband = 1 : length(Variable.nSubbands)
        nSubbands = Variable.nSubbands(iSubband);
        carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
        for iRealization = 1 : nRealizations
            channel = zeros(nTxs, nSubbands, nUsers);
            for iUser = 1 : nUsers
                % \boldsymbol{h}_{q,n}
                channel(:, :, iUser) = channel_tgn_e(pathloss, nSubbands, nTxs, carrierFrequency, fadingType);
            end
            % \boldsymbol{s_n}, v_{\text{out},q}
%             [waveformSu, voltageSu] = waveform_su(beta2, beta4, powerBudget, channel, tolerance);
            [waveformWsum, voltageWsum] = waveform_wsums(beta2, beta4, powerBudget, channel, tolerance, weight);
            [waveformCheWsum, voltageCheWsum] = waveform_che_wsum(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss);
            % % v_{\text{out},q}
            % voltageSu1 = harvester(beta2, beta4, waveformSu, channel);
            % voltageWsum1 = harvester(beta2, beta4, waveformWsum, channel);
        end
    end
end
voltageSu = mean(voltageSu, 3);
voltageWsum = mean(voltageWsum, 3);
save('data/wpt_wsum.mat');
%% Result
legendString = cell(2, length(Variable.nTxs));
figure('Name', sprintf('Average single user output voltage by SU WPT and WSum as a function of sinewaves'));
for iTx = 1 : length(Variable.nTxs)
    plot(Variable.nSubbands, voltageSu(iTx, :), 'Marker', x);
    legendString{1, iTx} = sprintf('SU WPT: M = %d', Variable.nTxs(iTx));
    plot(Variable.nSubbands, voltageWsum(iTx, :), 'Marker', o);
    legendString{2, iTx} = sprintf('WSum: M = %d', Variable.nTxs(iTx));
    hold on;
end
hold off;
grid minor;
legend(legendString);
xlabel('Number of tones')
ylabel('Average v_{out} [V]')
% savefig('results/wpt_wsum.fig');

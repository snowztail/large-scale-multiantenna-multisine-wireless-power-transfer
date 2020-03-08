function voltage = harvester(beta2, beta4, waveform, subchannel)


    [~, ~, nSubbands] = size(subchannel);
    subchannel = squeeze(subchannel);

    % the second order term in Taylor expansion
    term2 = 0;
    for iSubband = 1 : nSubbands
        term2 = term2 + waveform(:, iSubband)' * conj(subchannel(:, iSubband)) * subchannel(:, iSubband).' * waveform(:, iSubband);
    end
    % the fourth order term in Taylor expansion
    term4 = 0;
    for iSubband1 = 1 : nSubbands
        for iSubband2 = 1 : nSubbands
            for iSubband3 = 1 : nSubbands
                for iSubband4 = 1 : nSubbands
                    % output DC voltage if balanced
                    isBalanced = iSubband1 + iSubband2 == iSubband3 + iSubband4;
                    if isBalanced
                        term4 = term4 + (3 / 2) * waveform(:, iSubband3)' * conj(subchannel(:, iSubband3)) * subchannel(:, iSubband1).' * waveform(:, iSubband1) * waveform(:, iSubband4)' * conj(subchannel(:, iSubband4)) * subchannel(:, iSubband2).' * waveform(:, iSubband2);
                    end
                end
            end
        end
    end
    % truncate the voltage expression to the fourth order to capture fundamental behavior of rectifier nonlinearity
    voltage = real(beta2 * term2 + beta4 * term4);
end

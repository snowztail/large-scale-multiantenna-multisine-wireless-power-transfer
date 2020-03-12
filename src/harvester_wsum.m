function [sumVoltage, wsumVoltage] = harvester_wsum(beta2, beta4, waveform, channel, weight)
    % Function:
    %   - calculate the harvester output voltage
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - waveform [\boldsymbol{s_n}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - wsumVoltage [\sum w * v_{\text{out}}]: weighted sum of rectifier output DC voltage over all users
    %
    % Comment(s):
    %   - truncate the voltage expression to the fourth order to capture fundamental behavior of rectifier nonlinearity
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 08 Mar 20


    % single receive antenna
    [~, nSubbands, nUsers] = size(channel);

    % the second order term in Taylor expansion
    sumTerm2 = 0;
    wsumTerm2 = 0;
    for iSubband = 1 : nSubbands
        sumTerm2 = sumTerm2 + waveform(:, iSubband)' * conj(squeeze(channel(:, iSubband, :))) * squeeze(channel(:, iSubband, :)).' * waveform(:, iSubband);
        wsumTerm2 = wsumTerm2 + waveform(:, iSubband)' * conj(squeeze(channel(:, iSubband, :))) * diag(weight) * squeeze(channel(:, iSubband, :)).' * waveform(:, iSubband);
    end
    % the fourth order term in Taylor expansion
    sumTerm4 = 0;
    wsumTerm4 = 0;
    for iSubband1 = 1 : nSubbands
        for iSubband2 = 1 : nSubbands
            for iSubband3 = 1 : nSubbands
                for iSubband4 = 1 : nSubbands
                    % output DC voltage if balanced
                    isBalanced = iSubband1 + iSubband2 == iSubband3 + iSubband4;
                    if isBalanced
                        sumMiddleTerm = 0;
                        wsumMiddleTerm = 0;
                        for iUser = 1 : nUsers
                            sumMiddleTerm = sumMiddleTerm + conj(channel(:, iSubband3, iUser)) * channel(:, iSubband1, iUser).' * waveform(:, iSubband1) * waveform(:, iSubband4)' * conj(channel(:, iSubband4, iUser)) * channel(:, iSubband2, iUser).';
                            wsumMiddleTerm = wsumMiddleTerm + weight(iUser) * conj(channel(:, iSubband3, iUser)) * channel(:, iSubband1, iUser).' * waveform(:, iSubband1) * waveform(:, iSubband4)' * conj(channel(:, iSubband4, iUser)) * channel(:, iSubband2, iUser).';
                        end
                        sumTerm4 = sumTerm4 + (3 / 2) * waveform(:, iSubband3)' * sumMiddleTerm * waveform(:, iSubband2);
                        wsumTerm4 = wsumTerm4 + (3 / 2) * waveform(:, iSubband3)' * wsumMiddleTerm * waveform(:, iSubband2);
                    end
                end
            end
        end
    end
    sumVoltage = real(beta2 * sumTerm2 + beta4 * sumTerm4);
    wsumVoltage = real(beta2 * wsumTerm2 + beta4 * wsumTerm4);

end

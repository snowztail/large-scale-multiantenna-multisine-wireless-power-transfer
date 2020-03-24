function [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel)
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
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %   - minVoltage: minimum user voltage
    %
    % Comment(s):
    %   - truncate the voltage expression to the fourth order to capture fundamental behavior of rectifier nonlinearity
    %   - faster implementation based on an equivalent compact expression
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 14 Mar 20


    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    waveform = waveform(:);

    % \boldsymbol{M}_{q, k}
    channelMatrix = cell(nUsers, nSubbands);
    for iUser = 1 : nUsers
        subchannel = channel(:, :, iUser);
        % \boldsymbol{M}_{q}
        subchannelMatrix = conj(subchannel(:)) * subchannel(:).';
        for iSubband = 1 : nSubbands
            channelMatrix{iUser, iSubband} = zeros(nTxs * nSubbands);
            for jSubband = 1 : nSubbands + 1 - iSubband
                channelMatrix{iUser, iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = subchannelMatrix((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            end
        end
    end

    % v_{\text{out},q}
    userVoltage = zeros(1, nUsers);
    for iUser = 1 : nUsers
        userVoltage(iUser) = beta2 * waveform' * channelMatrix{iUser, 1} * waveform + (3 / 2) * beta4 * waveform' * channelMatrix{iUser, 1} * waveform * (waveform' * channelMatrix{iUser, 1} * waveform)';
        if nSubbands > 1
            for iSubband = 1 : nSubbands - 1
                userVoltage(iUser) = userVoltage(iUser) + 3 * beta4 * waveform' * channelMatrix{iUser, iSubband + 1} * waveform * (waveform' * channelMatrix{iUser, iSubband + 1} * waveform)';
            end
        end
    end
    userVoltage = real(userVoltage);
    % minimum user voltage
    minVoltage = min(userVoltage);
    % \sum v_{\text{out}}
    sumVoltage = sum(userVoltage);

end

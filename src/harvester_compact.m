function [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel)
    % Function:
    %   - calculate the harvester output voltage
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - truncate the voltage expression to the fourth order to capture fundamental behavior of rectifier nonlinearity
    %   - faster implementation based on an equivalent compact expression
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 14 Mar 20



    % * get block diagonal channel matrix
    [~, nSubbands, nUsers] = size(channel);
    % \boldsymbol{M}
    [matrixChannel] = matrix_channel(channel);

    % * compute output voltages
    waveformVector = vec(waveform);
    % v_{\text{out}}
    userVoltage = zeros(1, nUsers);
    for iUser = 1 : nUsers
        userVoltage(iUser) = real(beta2 * waveformVector' * matrixChannel{iUser, 1} * waveformVector + (3 / 2) * beta4 * waveformVector' * matrixChannel{iUser, 1} * waveformVector * (waveformVector' * matrixChannel{iUser, 1} * waveformVector)');
        if nSubbands > 1
            for iSubband = 1 : nSubbands - 1
                userVoltage(iUser) = userVoltage(iUser) + real(3 * beta4 * waveformVector' * matrixChannel{iUser, iSubband + 1} * waveformVector * (waveformVector' * matrixChannel{iUser, iSubband + 1} * waveformVector)');
            end
        end
    end
    % \min v_{\text{out}}
    minVoltage = min(userVoltage);
    % \sum v_{\text{out}}
    sumVoltage = sum(userVoltage);

end

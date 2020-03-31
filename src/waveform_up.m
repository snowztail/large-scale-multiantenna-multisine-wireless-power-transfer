function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_up(beta2, beta4, txPower, channel)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - for single-user and multi-user MISO systems
    %   - allocate power uniformly over all subbands
    %   - asymptotically optimal spatial beamformer
    %
    % Reference(s):
    %   - B. Clerckx and E. Bayguzina, "Waveform Design for Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 64, no. 23, pp. 6313–6328, Jan. 2016.
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20



    % * initialize complex carrier weight by uniform power allocation and spatial precoder by matched filter
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = repmat(sqrt(txPower / nSubbands / nUsers), [nSubbands, nUsers]);
    % \boldsymbol{w}
    precoder = zeros(nTxs, nSubbands, nUsers);
    for iUser = 1 : nUsers
        precoder(:, :, iUser) = conj(channel(:, :, iUser)) ./ vecnorm(channel(:, :, iUser), 2, 1);
    end

    % * construct waveform
    % \boldsymbol{s}
    waveform = sum(repmat(reshape(carrierWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* precoder, 3);
    % normalize waveform power
    waveform = sqrt(txPower) * waveform / norm(waveform, 'fro');

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

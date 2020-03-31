function [waveform, voltage] = waveform_ass(beta2, beta4, txPower, channel)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [h_n] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - voltage [v_{\text{out}}]: rectifier output DC voltage
    %
    % Comment(s):
    %   - for single-user MISO systems
    %   - allocate all power to the strongest subband
    %   - optimal for linear harvester model
    %
    % Reference(s):
    %   - B. Clerckx and E. Bayguzina, "Waveform Design for Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 64, no. 23, pp. 6313–6328, Jan. 2016.
    %
    % Author & Date: Yang (i@snowztail.com) - 08 Mar 20



    % * allocate all power to the strongest carrier
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * (max(vecnorm(channel, 2, 1)) == vecnorm(channel, 2, 1))';

    % * optimum single-user precoder is MRT
    % \boldsymbol{\tilde{s}}
    precoder = conj(channel) ./ vecnorm(channel, 2, 1);

    % * construct waveform
    % \boldsymbol{s_n}
    waveform = sum(repmat(reshape(carrierWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* precoder, 3);

    % * compute output voltage
    % v_{\text{out}}
    [voltage] = harvester_compact(beta2, beta4, waveform, channel);

end

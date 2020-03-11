function [waveform] = waveform_up(powerBudget, channel)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - powerBudget [P]: transmit power constraint
    %   - channel [h_{q, n}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s_n}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %
    % Comment(s):
    %   - for single-user MISO systems
    %   - allocate power uniformly over all subbands
    %
    % Reference(s):
    %   - B. Clerckx and E. Bayguzina, "Waveform Design for Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 64, no. 23, pp. 6313–6328, Jan. 2016.
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20

    % single receive antenna
    [~, nSubbands] = size(channel);
    % \boldsymbol{p}
    frequencyWeight = sqrt(powerBudget / nSubbands);
    % \boldsymbol{\tilde{s}_n}
    spatialPrecoder = conj(channel) ./ vecnorm(channel);
    % \boldsymbol{s_n}
    waveform = frequencyWeight * spatialPrecoder;

end

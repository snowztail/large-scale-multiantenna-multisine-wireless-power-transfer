function [waveform] = waveform_ass(powerBudget, subchannel)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - powerBudget [P]: transmit power constraint
    %   - subchannel [h_{q, n}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s_n}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
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


    % \boldsymbol{p}
    frequencyWeight = sqrt(powerBudget) * (max(vecnorm(subchannel)) == vecnorm(subchannel))';
    % \boldsymbol{\tilde{s}_n}
    spatialPrecoder = conj(subchannel) ./ vecnorm(subchannel);
    % \boldsymbol{s_n}
    waveform = frequencyWeight.' .* spatialPrecoder;

end

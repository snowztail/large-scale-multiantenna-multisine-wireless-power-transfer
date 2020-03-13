function [waveform, voltage] = waveform_up(beta2, beta4, powerBudget, channel, weight)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - powerBudget [P]: transmit power constraint
    %   - channel [h_{q, n}] (nTxs * nSubbands): channel frequency response at each subband
    %   - weight [w_q] (1 * nUsers): user weights
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %
    % Comment(s):
    %   - for single-user and multi-user MISO systems
    %   - allocate power uniformly over all subbands
    %
    % Reference(s):
    %   - B. Clerckx and E. Bayguzina, "Waveform Design for Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 64, no. 23, pp. 6313–6328, Jan. 2016.
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20


    % \boldsymbol{w}_n
    spatialPrecoder = sum(conj(channel) ./ vecnorm(channel, 2, 1), 3);
    % \boldsymbol{p}
    frequencyWeight = sqrt(powerBudget / norm(spatialPrecoder, 'fro') ^ 2);
    % \boldsymbol{s}_n
    waveform = frequencyWeight * spatialPrecoder;
    % v_{\text{out}, q}
    voltage = harvester(beta2, beta4, waveform, channel, weight);

end

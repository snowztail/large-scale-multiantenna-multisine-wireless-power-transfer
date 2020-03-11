function [channel] = channel_uncorrelated(pathloss, nSubbands, nTxs, fadingType)
    % Function:
    %   - simulate channel using the power delay profile of the IEEE TGn NLOS channel model E
    %
    % InputArg(s):
    %   - pathloss [\Lambda]: large-scale channel strength reduction
    %   - nSubbands [N]: number of subbands/subcarriers
    %   - nTxs [M]: number of transmit antennas
    %   - fadingType: "flat" or "selective"
    %
    % OutputArg(s):
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % Comment(s):
    %   - assume single receive antenna
    %   - the channel of a given user is sufficiently frequency-selective so that it is i.i.d. in space and frequency
    %   - the channel across users are fully uncorrelated
    %   - the reference pass loss is set to 60.046 dB for a distance of 10 m
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20

    nClusters = 4;
    nTaps = 18;
    tapPower = zeros(nTaps, nClusters);
    tapPower(:, 1) = db2pow([-2.6 -3.0 -3.5 -3.9 -4.5 -5.6 -6.9 -8.2 -9.8 -11.7 -13.9 -16.1 -18.3 -20.5 -22.9 -inf -inf -inf]');
    tapPower(:, 2) = db2pow([-inf -inf -inf -inf -1.8 -3.2 -4.5 -5.8 -7.1 -9.9 -10.3 -14.3 -14.7 -18.7 -19.9 -22.4 -inf -inf]');
    tapPower(:, 3) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -7.9 -9.6 -14.2 -13.8 -18.6 -18.1 -22.8 -inf -inf -inf]');
    tapPower(:, 4) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -20.6 -20.5 -20.7 -24.6]');

    % model taps as i.i.d. CSCG variables
    tapGain = repmat(sqrt(tapPower / 2), [1 1 nTxs]) .* (randn(nTaps, nClusters, nTxs) + 1i * randn(nTaps, nClusters, nTxs));
    tapGain = squeeze(sum(tapGain, 2));

    fading = zeros(nTxs, nSubbands);
    switch fadingType
    case 'selective'
        for iTx = 1 : nTxs
            for iSubband = 1 : nSubbands
                fading(iTx, iSubband) = sum(tapGain(:, iTx) .* exp(1i * 2 * pi * rand));
            end
        end
    case 'flat'
        for iTx = 1 : nTxs
            fading(iTx, :) = repmat(sum(tapGain(:, iTx) .* exp(1i * 2 * pi * rand)), [1 nSubbands]);
        end
    end
    channel = fading / sqrt(pathloss);

end

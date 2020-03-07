function [subchannel] = channel_tgn_e(distance, nSubbands, nRxs, nTxs, sampleFrequency, fadingType)


    nClusters = 4;
    nTaps = 18;
    tapDelay = 1e-9 * [0 10 20 30 50 80 110 140 180 230 280 330 380 430 490 560 640 730]';
    tapPower = zeros(nTaps, nClusters);
    tapPower(:, 1) = db2pow([-2.6 -3.0 -3.5 -3.9 -4.5 -5.6 -6.9 -8.2 -9.8 -11.7 -13.9 -16.1 -18.3 -20.5 -22.9 -inf -inf -inf]');
    tapPower(:, 2) = db2pow([-inf -inf -inf -inf -1.8 -3.2 -4.5 -5.8 -7.1 -9.9 -10.3 -14.3 -14.7 -18.7 -19.9 -22.4 -inf -inf]');
    tapPower(:, 3) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -7.9 -9.6 -14.2 -13.8 -18.6 -18.1 -22.8 -inf -inf -inf]');
    tapPower(:, 4) = db2pow([-inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -inf -20.6 -20.5 -20.7 -24.6]');

    % model taps as i.i.d. CSCG variables
    tapGain = cell(1, nClusters);
    for iCluster = 1 : nClusters
        tapGain{iCluster} = repmat(sqrt(tapPower(:, iCluster) / 2), [1 nRxs nTxs]) .* (randn(nTaps, nRxs, nTxs) + 1i * randn(nTaps, nRxs, nTxs));
    end
    tapGain = sum(cat(4, tapGain{:}), 4);

    fading = zeros(nSubbands, nRxs, nTxs);
    switch fadingType
    case "selective"
        for iRx = 1 : nRxs
            for iTx = 1 : nTxs
                for iSubband = 1 : nSubbands
                    fading(iSubband, iRx, iTx) = sum(tapGain(:, iRx, iTx) .* exp(1i * 2 * pi * sampleFrequency(iSubband) * tapDelay));
                end
            end
        end
    case "flat"
        for iRx = 1 : nRxs
            for iTx = 1 : nTxs
                fading(:, iRx, iTx) = repmat(sum(tapGain(:, iRx, iTx) .* exp(1i * 2 * pi * mean(sampleFrequency) * tapDelay)), [nSubbands 1]);
            end
        end
    end

    % large-scale fading
    pathlossExponent = 2;
    % shadowingStdDevDb = 3;
    pathloss = db2pow(60.046 + 10 * pathlossExponent * log10(distance / 10));
    % shadowing = db2pow(shadowingStdDevDb * randn);

    % subchannel = fading / pathloss / shadowing;
    subchannel = fading / pathloss;

end

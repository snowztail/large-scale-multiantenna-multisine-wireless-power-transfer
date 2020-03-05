function frequencyPrecoder = wpt_su(beta4, powerBudget, nSubbands, nTxs, subchannel)

    % \tilde{s}
    spatialPrecoder = zeros(nTxs * nSubbands, 1);
    % \tilde{s}_n
    spatialSubprecoder = reshape(spatialPrecoder, [nTxs, nSubbands]);

    % ? initialize p by uniform power allocation
    frequencyPrecoder = sqrt(powerBudget / nSubbands) * ones(nSubbands, 1);
    for iSubband = 1 : nSubbands
        spatialSubprecoder(:, iSubband) = conj(subchannel(:, iSubband)) / norm(subchannel(:, iSubband));
    end

    % \boldsymbol{X}
    frequencyPrecoderMatrix = frequencyPrecoder * frequencyPrecoder';

    % \boldsymbol{M}''
    channelNormMatrix = vecnorm(subchannel).' * vecnorm(subchannel);
    % \boldsymbol{M}''_k
    channelNormDiagMatrix = cell(nSubbands, 1);
    % \boldsymbol{t}
    auxiliary = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        channelNormDiagMatrix{iSubband} = diag(diag(channelNormMatrix, iSubband - 1), iSubband - 1);
        auxiliary(iSubband) = trace(channelNormDiagMatrix{iSubband} * frequencyPrecoderMatrix);
    end
    % \boldsymbol{A}_0
    coefficient = diag(beta4 * [-1.5; -3 * ones(nSubbands - 1, 1)]);
    % g_{boldsymbol{t}}
    target4 = auxiliary' * coefficient * auxiliary;
end

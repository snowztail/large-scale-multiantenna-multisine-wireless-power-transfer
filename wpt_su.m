function [weight] = wpt_su(beta2, beta4, powerBudget, nSubbands, nTxs, subchannel, threshold)



    %% Initialization
    % ? initialize p by uniform power allocation
    frequencyWeight = sqrt(powerBudget / nSubbands) * ones(nSubbands, 1);
    % \boldsymbol{X}
    frequencyWeightMatrix = frequencyWeight * frequencyWeight';
    % \boldsymbol{M}''
    channelNormMatrix = vecnorm(subchannel).' * vecnorm(subchannel);
    % \boldsymbol{M}''_k
    channelNormDiagMatrix = cell(nSubbands, 1);
    % \boldsymbol{t}
    auxiliary = zeros(nSubbands, 1);
    for iSubband = 1 : nSubbands
        channelNormDiagMatrix{iSubband} = diag(diag(channelNormMatrix, iSubband - 1), iSubband - 1);
        auxiliary(iSubband) = trace(channelNormDiagMatrix{iSubband} * frequencyWeightMatrix);
    end
    %% Iteration
    isConverged = false;
    while ~isConverged
        % \boldsymbol{C}''_1
        termC1 = - (beta2 + 3 * beta4 * auxiliary(1)) / 2 * channelNormDiagMatrix{1} - 3 * beta4 * sum(cat(3, channelNormDiagMatrix{2 : end}) .* reshape(conj(auxiliary(2 : end)), [1, 1, nSubbands - 1]), 3);
        % \boldsymbol{A}''_1
        termA1 = termC1 + termC1';

        % * Solve rank-1 \boldsymbol{X}^{\star} in closed form (low complexity)
        % \boldsymbol{x}^{\star}
        [v, d] = eig(termA1);
        frequencyWeight = sqrt(powerBudget) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        frequencyPrecoderMatrix_ = frequencyWeight * frequencyWeight';
        % Update \boldsymbol{t}
        auxiliary = zeros(nSubbands, 1);
        for iSubband = 1 : nSubbands
            auxiliary(iSubband) = trace(channelNormDiagMatrix{iSubband} * frequencyPrecoderMatrix_);
        end
        % test convergence
        if (norm(frequencyPrecoderMatrix_, 'fro') - norm(frequencyWeightMatrix, 'fro')) / norm(frequencyPrecoderMatrix_, 'fro') <= threshold
            isConverged = true;
        end
    end
    spatialPrecoder = conj(subchannel) ./ vecnorm(subchannel);
    % s_n
    weight = zeros(nTxs, nSubbands);
    for iSubband = 1 : nSubbands
        weight(:, iSubband) = frequencyWeight(iSubband) * spatialPrecoder(:, iSubband);
    end

end

% % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
% cvx_begin
%     cvx_solver mosek
%     variable frequencyWeightMatrix
%     minimize trace(termA1 * frequencyWeightMatrix)
%     subject to
%         trace(frequencyWeightMatrix) <= powerBudget;
% cvx_end
% % update auxiliary terms
% auxiliary = zeros(nSubbands, 1);
% for iSubband = 1 : nSubbands
%     auxiliary(iSubband) = trace(channelNormDiagMatrix{iSubband} * frequencyWeightMatrix);
% end

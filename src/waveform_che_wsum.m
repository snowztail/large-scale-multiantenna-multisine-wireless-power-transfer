function [waveform, asymWaveform] = waveform_che_wsum(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss)

    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    % ? users have the same pathloss
    % pathloss = rand(1, nUsers);
    pathloss = pathloss * ones(1, nUsers);
    % ? initialize \boldsymbol{p}_q by uniform power allocation over subbands of a given user (the power across users depend on pathloss)
    frequencyWeight = sqrt(ones(nSubbands, nUsers) / nSubbands / nUsers ./ pathloss);

    % \boldsymbol{M}'_{k}
    shiftMatrix = cell(1, nSubbands);
    for iSubband = 1 : nSubbands
        shiftMatrix{iSubband} = diag(diag(ones(nSubbands), iSubband - 1), iSubband - 1);
    end
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        for iSubband = 1 : nSubbands
            auxiliary(iUser, iSubband) = frequencyWeight(:, iUser)' * shiftMatrix{iSubband} * frequencyWeight(:, iUser);
        end
    end

    isConverged = false;
    while ~isConverged
        % \boldsymbol{C}'_{q, 1}
        termC1 = cell(nUsers, 1);
        % \boldsymbol{A}'_{q, 1}
        termA1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            termC1{iUser} = - (beta2 * powerBudget * nTxs * pathloss(iUser) ^ 2 + 3 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * shiftMatrix{1} - 3 * beta4 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, shiftMatrix{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            termA1{iUser} = weight(iUser) * (termC1{iUser} + termC1{iUser}');
        end
        % \boldsymbol{A}'_1
        termA1 = blkdiag(termA1{:});
        % \boldsymbol{\Lambda}
        pathlossMatrix = diag(repelem(pathloss, nSubbands));

        % * Solve \bar{\boldsymbol{p}} in closed form (low complexity)
        [v, d] = eig(termA1 / pathlossMatrix);
        vMin = v(:, diag(d) == min(diag(d)));
        % \bar{\boldsymbol{p}}
        frequencyWeight_ = sqrt(1 / (vMin' * pathlossMatrix * vMin)) * vMin;
        clearvars v d vMin;
        % reshape for visualization
        frequencyWeight_ = reshape(frequencyWeight_, nSubbands, nUsers);

        % update \boldsymbol{t}_q
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = frequencyWeight_(:, iUser)' * shiftMatrix{iSubband} * frequencyWeight_(:, iUser);
            end
        end
        % test convergence
        if (norm(frequencyWeight_ - frequencyWeight, 'fro')) / norm(frequencyWeight_, 'fro') <= tolerance
            isConverged = true;
        end
        frequencyWeight = frequencyWeight_;
    end
    % \bar{\boldsymbol{s}}_n
    waveform = sum(repmat(reshape(frequencyWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* conj(channel), 3) / sqrt(nTxs);
    % \boldsymbol{s}_{\text{asym}}
    asymWaveform = sqrt(powerBudget) * waveform / norm(waveform, 'fro');

end

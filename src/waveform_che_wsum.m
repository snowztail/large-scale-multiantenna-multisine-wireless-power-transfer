function [waveform, sumVoltage, userVoltage] = waveform_che_wsum(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - powerBudget [P]: transmit power constraint
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w_q] (1 * nUsers): user weights
    %   - pathloss [\Lambda] (1 * nUsers): user pathlosses
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_{\text{asym}}] (nTxs * nSubbands): the asymptotically optimal complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %
    % Comment(s):
    %   - for multi-user MISO systems
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - note the notation difference with the single-user case
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 12 Mar 20


    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    % ? users have the same pathloss
    % pathloss = rand(1, nUsers);
    pathloss = 1 / pathloss * ones(1, nUsers);
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
            termC1{iUser} = - ((beta2 * powerBudget * nTxs * pathloss(iUser) ^ 2 + 3 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * shiftMatrix{1});
            if nSubbands > 1
                termC1{iUser} = termC1{iUser} - 3 * beta4 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, shiftMatrix{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            termA1{iUser} = weight(iUser) * (termC1{iUser} + termC1{iUser}');
        end
        % \boldsymbol{A}'_1
        termA1 = blkdiag(termA1{:});
        % \boldsymbol{\Lambda}
        pathlossMatrix = diag(repelem(pathloss, nSubbands));

        % * Solve \bar{\boldsymbol{p}} in closed form (low complexity)
        [v, d] = eig(pathlossMatrix \ termA1);
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

    % % * asymptotic out voltage v_{\text{out},q}'
    % voltage = zeros(1, nUsers);
    % for iUser = 1 : nUsers
    %     voltage(iUser) = beta2 * powerBudget * nTxs * pathloss(iUser) ^ 2 * frequencyWeight(:, iUser)' * frequencyWeight(:, iUser) + 3 / 2 * beta4 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * (frequencyWeight(:, iUser)' * shiftMatrix{1} * frequencyWeight(:, iUser)) * (frequencyWeight(:, iUser)' * shiftMatrix{1} * frequencyWeight(:, iUser))';
    %     if nSubbands > 1
    %         for iSubband = 1 : nSubbands - 1
    %             voltage(iUser) = voltage(iUser) + 3 * beta4 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * (frequencyWeight(:, iUser)' * shiftMatrix{iSubband + 1} * frequencyWeight(:, iUser)) * (frequencyWeight(:, iUser)' * shiftMatrix{iSubband + 1} * frequencyWeight(:, iUser))';
    %         end
    %     end
    % end
    % % \sum v_{\text{out}}'
    % sumVoltage = real(sum(voltage));
    % % \sum w * v_{\text{out}}'
    % wsumVoltage = real(sum(weight .* voltage));

    % \bar{\boldsymbol{s}}_n
    normalizedWaveform = sum(repmat(reshape(frequencyWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* conj(channel), 3) / sqrt(nTxs);
    % \boldsymbol{s}_{\text{asym}}
    waveform = sqrt(powerBudget) * normalizedWaveform / norm(normalizedWaveform, 'fro');
    % \sum v_{\text{out}}, v\{\text{out}, q}
    [sumVoltage, userVoltage] = harvester_compact(beta2, beta4, waveform, channel);
end

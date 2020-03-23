function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_che_rr(beta2, beta4, powerBudget, channel, tolerance, weight, pathloss)
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
    %   - minVoltage: minimum user voltage
    %
    % Comment(s):
    %   - maximize the minimum user voltage
    %   - for multi-user MISO systems with number of user no larger than 3
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - in each iteration, we first obtain the high rank covariance matrix by CVX, then perform rank reduction for a rank-1 solution
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 17 Mar 20


    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    % ? users have the same pathloss
    % pathloss = rand(1, nUsers);
    pathloss = 1 ./ pathloss;
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
    % \boldsymbol{A}_0
    termA0 = diag(-3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        % \bar{c}_q'
        termBarC = zeros(nUsers, 1);
        % \boldsymbol{C}'_{q, 1}
        termC1 = cell(nUsers, 1);
        % \boldsymbol{A}'_{q, 1}
        termA1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            termBarC(iUser) = - real(conj(auxiliary(iUser, :)) * termA0 * auxiliary(iUser, :).' * nTxs ^ 2 * powerBudget ^ 2 * pathloss(iUser) ^ 4);
            termC1{iUser} = - ((beta2 * powerBudget * nTxs * pathloss(iUser) ^ 2 + 3 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * shiftMatrix{1});
            if nSubbands > 1
                termC1{iUser} = termC1{iUser} - 3 * beta4 * powerBudget ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, shiftMatrix{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            termA1{iUser} = weight(iUser) * (termC1{iUser} + termC1{iUser}');
        end

        % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
        cvx_begin quiet
            % \boldsymbol{X}
            variable highRankFrequencyWeightMatrix(nSubbands, nSubbands, nUsers) complex semidefinite;
            target = cvx(zeros(1, nUsers));
            sumPower = 0;
            for iUser = 1 : nUsers
                target(iUser) = trace(termA1{iUser} * highRankFrequencyWeightMatrix(:, :, iUser)) + termBarC(iUser);
                sumPower = sumPower + pathloss(iUser) * trace(highRankFrequencyWeightMatrix(:, :, iUser));
            end
            minimize(max(target));
            subject to
                sumPower = 1;
        cvx_end
        [~, userIndex] = max(target);
    end

    minVoltage = 1;
end
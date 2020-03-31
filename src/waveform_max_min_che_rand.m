function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_che_rand(beta2, beta4, txPower, channel, tolerance, weight, pathloss)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w_q] (1 * nUsers): user weights
    %   - pathloss [\Lambda] (1 * nUsers): user pathlosses
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_{\text{asym}}] (nTxs * nSubbands): the asymptotically optimal complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - maximize the minimum user voltage
    %   - for multi-user MISO systems with arbitrary number of user
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - we first obtain the high rank covariance matrix by CVX, then generate waveform candidates based on random vectors and pick the optimal one as a rank-1 solution
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 25 Mar 20



    [nTxs, nSubbands, nUsers] = size(channel);
    % ? users have the same pathloss
    % pathloss = rand(1, nUsers);
    pathloss = 1 ./ pathloss;
    % % ? initialize \boldsymbol{p}_q by uniform power allocation over subbands of a given user (the power across users depend on pathloss)
    % carrierWeight = sqrt(ones(nSubbands, nUsers) / nSubbands / nUsers ./ pathloss);

    carrierWeight = squeeze(sum(channel)) / norm(squeeze(sum(channel)), 'fro') * sqrt(txPower);

    % \boldsymbol{M}'_{k}
    matrixShift = cell(1, nSubbands);
    for iSubband = 1 : nSubbands
        matrixShift{iSubband} = diag(diag(ones(nSubbands), iSubband - 1), iSubband - 1);
    end
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    % \boldsymbol{X}_{q}
    carrierWeightMatrix = zeros(nSubbands, nSubbands, nUsers);
    for iUser = 1 : nUsers
        carrierWeightMatrix(:, :, iUser) = carrierWeight(:, iUser) * carrierWeight(:, iUser)';
        for iSubband = 1 : nSubbands
            auxiliary(iUser, iSubband) = trace(matrixShift{iSubband} * carrierWeightMatrix(:, :, iUser));
        end
    end

    isConverged = false;
    counter = 0;
    maxTarget = 0;
    % \boldsymbol{A}_0
    termA0 = diag(- 3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        counter = counter + 1;
        % \bar{c}_q'
        termBarC = zeros(nUsers, 1);
        % \boldsymbol{C}'_{q, 1}
        C1 = cell(nUsers, 1);
        % \boldsymbol{A}'_{q, 1}
        A1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            termBarC(iUser) = - real(conj(auxiliary(iUser, :)) * termA0 * auxiliary(iUser, :).' * nTxs ^ 2 * txPower ^ 2 * pathloss(iUser) ^ 4);
            C1{iUser} = - ((beta2 * txPower * nTxs * pathloss(iUser) ^ 2 + 3 * txPower ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * matrixShift{1});
            if nSubbands > 1
                C1{iUser} = C1{iUser} - 3 * beta4 * txPower ^ 2 * nTxs ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, matrixShift{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end

        % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
        cvx_begin quiet
            % \boldsymbol{X}
            variable highRankcarrierWeightMatrix(nSubbands, nSubbands, nUsers) complex semidefinite;
            target = cvx(zeros(1, nUsers));
            traceSum = 0;
            for iUser = 1 : nUsers
                target(iUser) = trace(A1{iUser} * highRankcarrierWeightMatrix(:, :, iUser)) + termBarC(iUser);
                traceSum = traceSum + pathloss(iUser) * trace(highRankcarrierWeightMatrix(:, :, iUser));
            end
            minimize(max(target));
            subject to
                traceSum == 1;
        cvx_end
        [~, userIndex] = max(target);

        % * Derive a best rank-1 solution from randomly generated vectors
        carrierWeightMatrix_ = highRankcarrierWeightMatrix;
        % denote term \boldsymbol{A}'_{q, 1} as \boldsymbol{B}_{1, q} for any q ~= q_0
        B1 = A1;
        B1{userIndex} = 0;
%         % denote term - \boldsymbol{A}'_{q_0, 1} as \boldsymbol{B}_{1, q_0}
%         B1Max = - A1{userIndex};
        % denote term \Lambda_q * \boldsymbol{I} as \boldsymbol{B}_{2, q}
        B2 = cell(nUsers, 1);
        termQ = cell(nUsers, 1);
        carrierWeightSqrtMatrix = cell(nUsers, 1);
        for iUser = 1 : nUsers
            B2{iUser} = pathloss(iUser) * eye(nSubbands);
            carrierWeightSqrtMatrix{iUser} = carrierWeightMatrix_(:, :, iUser) ^ (1 / 2);
            [v, ~] = eig(carrierWeightSqrtMatrix{iUser} * B1{iUser} * carrierWeightSqrtMatrix{iUser});
            termQ{iUser} = v' * carrierWeightSqrtMatrix{iUser} * B2{iUser} * carrierWeightSqrtMatrix{iUser} * v;
            randomVector = randomize(termQ{iUser});
            carrierWeight(:, iUser) = carrierWeightSqrtMatrix{iUser} * v * randomVector;
            carrierWeightMatrix_(:, :, iUser) = carrierWeight(:, iUser) * carrierWeight(:, iUser)';
            clearvars v;
        end

        % Update \boldsymbol{t}_{q, k}
        auxiliary = zeros(nUsers, nSubbands);
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixShift{iSubband} * carrierWeightMatrix_(:, :, iUser));
            end
        end

        % Update target function
        target = zeros(nUsers, 1);
        for iUser = 1 : nUsers
            target(iUser) = real(trace(A1{iUser} * carrierWeightMatrix_(:, :, iUser)) + termBarC(iUser));
        end
        maxTarget_ = max(target);

        % test convergence
        temp = abs(maxTarget_ - maxTarget) / abs(maxTarget_)
        if abs(maxTarget_ - maxTarget) / abs(maxTarget_) <= tolerance || counter >= 1e2
            isConverged = true;
        end
        maxTarget = maxTarget_;

    end
    % \bar{\boldsymbol{s}}_n
    normalizedWaveform = sum(repmat(reshape(carrierWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* conj(channel), 3) / sqrt(nTxs);
    % \boldsymbol{s}_{\text{asym}}
    waveform = sqrt(txPower) * normalizedWaveform / norm(normalizedWaveform, 'fro');
    % \sum v_{\text{out}}, v\{\text{out}, q}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_rand(beta2, beta4, txPower, channel, tolerance, weight, nCandidates)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the minimum voltage with randomized power allocation
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w_q] (1 * nUsers): user weights
    %   - nCandidates: number of random feasible rank-1 solutions to generate
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - maximize the minimum user voltage
    %   - for multi-user MISO systems with arbitrary number of user
    %   - we first obtain the high rank covariance matrix by CVX, then generate waveform candidates based on random vectors and pick the optimal one as a rank-1 solution
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 17 Mar 20



    [nTxs, nSubbands, nUsers] = size(channel);
    % ? initialize \boldsymbol{s} by uniform power allocation and omidirectional beamforming
    waveform = sqrt(txPower / nTxs / nSubbands) * ones(nTxs, nSubbands);
    % \boldsymbol{X}
    waveformMatrix = waveform(:) * waveform(:)';

    % \boldsymbol{M}_{q, k}
    matrixChannel = cell(nUsers, nSubbands);
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        subchannel = channel(:, :, iUser);
        % \boldsymbol{M}_{q}
        submatrixChannel = conj(subchannel(:)) * subchannel(:).';
        for iSubband = 1 : nSubbands
            matrixChannel{iUser, iSubband} = zeros(nTxs * nSubbands);
            for jSubband = 1 : nSubbands + 1 - iSubband
                matrixChannel{iUser, iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = submatrixChannel((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            end
            auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * waveformMatrix);
        end
    end

    isConverged = false;
    % \boldsymbol{A}_0
    A0 = diag(-3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        % \bar{c}_q
        cBar = zeros(1, nUsers);
        % \boldsymbol{C}_{q, 1}
        C1 = cell(1, nUsers);
        % \boldsymbol{A}_{q, 1}
        A1 = cell(1, nUsers);
        for iUser = 1 : nUsers
            cBar(iUser) = - real(conj(auxiliary(iUser, :)) * A0 * auxiliary(iUser, :).');
            C1{iUser} = - (beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * matrixChannel{iUser, 1};
            if nSubbands > 1
                C1{iUser} = C1{iUser} - weight(iUser) * 3 * beta4 * sum(cat(3, matrixChannel{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end

        % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
        cvx_begin quiet
            % \boldsymbol{X}
            variable highRankWaveformMatrix(nTxs * nSubbands, nTxs * nSubbands) complex semidefinite;
            target = cvx(zeros(1, nUsers));
            for iUser = 1 : nUsers
                target(iUser) = trace(A1{iUser} * highRankWaveformMatrix) + cBar(iUser);
            end
            minimize(max(target));
            subject to
                trace(highRankWaveformMatrix) <= txPower;
        cvx_end
        waveformMatrix_ = highRankWaveformMatrix;

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * waveformMatrix_);
            end
        end

        % test convergence
        if (norm(waveformMatrix_ - waveformMatrix, 'fro')) / norm(waveformMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        waveformMatrix = waveformMatrix_;
    end

    % * Derive a best rank-1 solution from randomly generated vectors
    [v, d] = eig(waveformMatrix);
    % generate waveform candidates based on random vectors
    candidateWaveform = cell(nCandidates, 1);
    candidateTarget = zeros(nCandidates, nUsers);
    for iCandidate = 1 : nCandidates
        % candidateWaveform{iCandidate} = v * d .^ (1 / 2) * (sqrt(rand(nTxs * nSubbands, 1)) .* (exp(1i * 2 * pi * rand(nTxs * nSubbands, 1))));
        candidateWaveform{iCandidate} = v * d .^ (1 / 2) * (exp(1i * 2 * pi * rand(nTxs * nSubbands, 1)));
        candidateWaveformMatrix = candidateWaveform{iCandidate} * candidateWaveform{iCandidate}';

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * candidateWaveformMatrix);
            end
        end

        % \bar{c}_q
        cBar = zeros(1, nUsers);
        % \boldsymbol{C}_{q, 1}
        C1 = cell(1, nUsers);
        % \boldsymbol{A}_{q, 1}
        A1 = cell(1, nUsers);
        for iUser = 1 : nUsers
            cBar(iUser) = - real(conj(auxiliary(iUser, :)) * A0 * auxiliary(iUser, :).');
            C1{iUser} = - (beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * matrixChannel{iUser, 1};
            if nSubbands > 1
                C1{iUser} = C1{iUser} - weight(iUser) * 3 * beta4 * sum(cat(3, matrixChannel{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end

        for iUser = 1 : nUsers
            candidateTarget(iCandidate, iUser) = real(trace(A1{iUser} * candidateWaveformMatrix) + cBar(iUser));
        end

    end
    [~, candidateIndex] = min(max(candidateTarget, [], 2));

    waveform = candidateWaveform{candidateIndex};

    % v_{\text{out}, q}
    userVoltage = zeros(1, nUsers);
    for iUser = 1 : nUsers
        userVoltage(iUser) = beta2 * waveform' * matrixChannel{iUser, 1} * waveform + (3 / 2) * beta4 * waveform' * matrixChannel{iUser, 1} * waveform * (waveform' * matrixChannel{iUser, 1} * waveform)';
        if nSubbands > 1
            for iSubband = 1 : nSubbands - 1
                userVoltage(iUser) = userVoltage(iUser) + 3 * beta4 * waveform' * matrixChannel{iUser, iSubband + 1} * waveform * (waveform' * matrixChannel{iUser, iSubband + 1} * waveform)';
            end
        end
    end
    userVoltage = real(userVoltage);
    minVoltage = min(userVoltage);
    % \boldsymbol{s_n}
    waveform = reshape(waveform, [nTxs, nSubbands]);
    % \sum v_{\text{out}}
    sumVoltage = sum(userVoltage);

end

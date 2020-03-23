function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_rand(beta2, beta4, powerBudget, channel, tolerance, weight, nCandidates)
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
    %   - nCandidates: number of random feasible rank-1 solutions to generate
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %   - minVoltage: minimum user voltage
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


    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    % ? initialize \boldsymbol{s} by uniform power allocation and omidirectional beamforming
    waveform = sqrt(powerBudget / nTxs / nSubbands) * ones(nTxs, nSubbands);
    % \boldsymbol{X}
    waveformMatrix = waveform(:) * waveform(:)';

    % \boldsymbol{M}_{q, k}
    channelMatrix = cell(nUsers, nSubbands);
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        subchannel = channel(:, :, iUser);
        % \boldsymbol{M}_{q}
        subchannelMatrix = conj(subchannel(:)) * subchannel(:).';
        for iSubband = 1 : nSubbands
            channelMatrix{iUser, iSubband} = zeros(nTxs * nSubbands);
            for jSubband = 1 : nSubbands + 1 - iSubband
                channelMatrix{iUser, iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = subchannelMatrix((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            end
            auxiliary(iUser, iSubband) = trace(channelMatrix{iUser, iSubband} * waveformMatrix);
        end
    end

    isConverged = false;
    % \boldsymbol{A}_0
    termA0 = diag(-3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        % \bar{c}_q
        termBarC = zeros(1, nUsers);
        % \boldsymbol{C}_{q, 1}
        termC1 = cell(1, nUsers);
        % \boldsymbol{A}_{q, 1}
        termA1 = cell(1, nUsers);
        for iUser = 1 : nUsers
            termBarC(iUser) = - real(conj(auxiliary(iUser, :)) * termA0 * auxiliary(iUser, :).');
            termC1{iUser} = - ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * channelMatrix{iUser, 1});
            if nSubbands > 1
                termC1{iUser} = termC1{iUser} - weight(iUser) * 3 * beta4 * sum(cat(3, channelMatrix{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            termA1{iUser} = termC1{iUser} + termC1{iUser}';
        end

        % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
        cvx_begin quiet
            % \boldsymbol{X}
            variable highRankWaveformMatrix(nTxs * nSubbands, nTxs * nSubbands) complex semidefinite;
            target = cvx(zeros(1, nUsers));
            for iUser = 1 : nUsers
                target(iUser) = trace(termA1{iUser} * highRankWaveformMatrix) + termBarC(iUser);
            end
            minimize(max(target));
            subject to
                trace(highRankWaveformMatrix) <= powerBudget;
        cvx_end
        waveformMatrix_ = highRankWaveformMatrix;

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(channelMatrix{iUser, iSubband} * waveformMatrix_);
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
        radius = rand(nTxs * nSubbands, 1);
        phase = 2 * pi * rand(nTxs * nSubbands, 1);
        candidateWaveform{iCandidate} = v * d ^ (1 / 2) * (sqrt(radius) .* (cos(phase) + 1i * sin(phase)));
        candidateWaveformMatrix = candidateWaveform{iCandidate} * candidateWaveform{iCandidate}';

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(channelMatrix{iUser, iSubband} * candidateWaveformMatrix);
            end
        end

        % \bar{c}_q
        termBarC = zeros(1, nUsers);
        % \boldsymbol{C}_{q, 1}
        termC1 = cell(1, nUsers);
        % \boldsymbol{A}_{q, 1}
        termA1 = cell(1, nUsers);
        for iUser = 1 : nUsers
            termBarC(iUser) = - real(conj(auxiliary(iUser, :)) * termA0 * auxiliary(iUser, :).');
            termC1{iUser} = - ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * channelMatrix{iUser, 1});
            if nSubbands > 1
                termC1{iUser} = termC1{iUser} - weight(iUser) * 3 * beta4 * sum(cat(3, channelMatrix{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            termA1{iUser} = termC1{iUser} + termC1{iUser}';
        end

        for iUser = 1 : nUsers
            candidateTarget(iCandidate, iUser) = real(trace(termA1{iUser} * candidateWaveformMatrix) + termBarC(iUser));
        end

    end
    [~, candidateIndex] = min(max(candidateTarget, [], 2));

    waveform = candidateWaveform{candidateIndex};

    % v_{\text{out}, q}
    userVoltage = zeros(1, nUsers);
    for iUser = 1 : nUsers
        userVoltage(iUser) = beta2 * waveform' * channelMatrix{iUser, 1} * waveform + (3 / 2) * beta4 * waveform' * channelMatrix{iUser, 1} * waveform * (waveform' * channelMatrix{iUser, 1} * waveform)';
        if nSubbands > 1
            for iSubband = 1 : nSubbands - 1
                userVoltage(iUser) = userVoltage(iUser) + 3 * beta4 * waveform' * channelMatrix{iUser, iSubband + 1} * waveform * (waveform' * channelMatrix{iUser, iSubband + 1} * waveform)';
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
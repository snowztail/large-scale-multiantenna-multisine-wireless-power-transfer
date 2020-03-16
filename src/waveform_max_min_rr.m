function [waveform, sumVoltage, userVoltage] = waveform_max_min_rr(beta2, beta4, powerBudget, channel, tolerance, weight)
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
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}, q}]: individual user voltages
    %
    % Comment(s):
    %   - for multi-user MISO systems
    %   - the frequency domain power allocation and spatial beamforming are coupled for multi-user scenario
    %   - the suboptimal spatial beamformer is obtained in a closed form based on a linear harvester model
    %   - note the notation difference with the single-user case
    %   - obtain the rank-1 frequency weight matrix in closed form for a complexity of \mathcal{O}((N) ^ 3)
    %   - as all eigenvalues of \boldsymbol{A}'''_1 are negative, we choose the eigenvector corresponding to the minimum eigenvalue (maximum in magnitude)
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 15 Mar 20

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
            variable waveformMatrix(nTxs * nSubbands, nTxs * nSubbands) complex semidefinite;
            target = cvx(zeros(1, nUsers));
            for iUser = 1 : nUsers
                target(iUser) = trace(termA1{iUser} * waveformMatrix) + termBarC(iUser);
            end
            minimize(max(target));
            subject to
                trace(waveformMatrix) <= powerBudget;
        cvx_end
        [~, userIndex] = max(target);

        % * Rank reduction for separable SDP
        waveformRank = rank(waveformMatrix);
        while waveformRank ^ 2 > nUsers
            % decompose waveform matrix as a product of a matrix V and its Hermitian
            waveformFactor = cholcov(waveformMatrix)';
            % flatten trace equations to standard linear equations
            coefficient = zeros(nUsers, size(waveformFactor, 2) ^ 2);
            for iUser = 1 : nUsers
                coefficient(iUser, :) = reshape((waveformFactor' * (termA1{iUser} - termA1{userIndex}) * waveformFactor).', 1, []);
            end
            % obtain an orthonormal basis for the null space of the coefficient matrix
            delta = null(coefficient);
            % nonzero solution can be obtained as a linear combination of null space basis vectors
            delta = reshape(delta(:, 1), [waveformRank, waveformRank]);
            delta = (delta + delta') / 2;
            [~, d] = eig(delta);
            dominantEigenvalue = d(abs(diag(d)) == max(abs(diag(d))));
            waveformMatrix = waveformFactor * (eye(waveformRank) - 1 / dominantEigenvalue(1) * delta) * waveformFactor';
            % fix accuracy issue and ensure strict symmetric
            waveformMatrix = (waveformMatrix + waveformMatrix') / 2;
            % update waveform rank
            waveformRank = rank(waveformMatrix);
        end




    end


    sumVoltage = 1;
    userVoltage = 1;
end

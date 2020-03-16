clear; close all; clc; initialize; config_max_min_rr;
channel = channel_tgn_e(pathloss, nTxs, nSubbands, nUsers, carrierFrequency, fadingType);
% [waveform, sumVoltage, userVoltage] = waveform_max_min_rr(beta2, beta4, powerBudget, channel, tolerance, weight);
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
            delta = reshape(delta(:, 1), [size(waveformFactor, 2), size(waveformFactor, 2)]);
            delta = (delta + delta') / 2;
            [~, d] = eig(delta);
            dominantEigenvalue = d(abs(diag(d)) == max(abs(diag(d))));
            waveformMatrix = waveformFactor * (eye(size(waveformFactor, 2)) - 1 / dominantEigenvalue(1) * delta) * waveformFactor';
            % fix accuracy issue and ensure strict symmetric
            waveformMatrix = (waveformMatrix + waveformMatrix') / 2;
            % update waveform rank
            waveformRank = rank(waveformMatrix);
        end
    end
    

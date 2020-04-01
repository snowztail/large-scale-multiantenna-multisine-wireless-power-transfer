function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_rr(beta2, beta4, txPower, channel, tolerance, weight)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the minimum voltage with rank reduction
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w] (1 * nUsers): user weights
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - maximize the minimum user voltage
    %   - for multi-user MISO systems with number of user no larger than 3
    %   - in each iteration, we first obtain the high rank covariance matrix by CVX, then perform rank reduction for a rank-1 solution
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %   - Y. Huang and D. Palomar, "Rank-Constrained Separable Semidefinite Programming With Applications to Optimal Beamforming," IEEE Transactions on Signal Processing, vol. 58, no. 2, pp. 664–678, 2010.
    %
    % Author & Date: Yang (i@snowztail.com) - 17 Mar 20



    % * initialize complex carrier weight by channel strength and spatial precoder by matched filter
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * permute(sum(channel, 1), [2, 3, 1]) / norm(squeeze(sum(channel)), 'fro');
    % \boldsymbol{\tilde{s}}
    precoder = zeros(nTxs, nSubbands, nUsers);
    for iUser = 1 : nUsers
        precoder(:, :, iUser) = conj(channel(:, :, iUser)) ./ vecnorm(channel(:, :, iUser), 2, 1);
    end
    % \boldsymbol{s}
    waveform = sum(permute(carrierWeight, [3, 1, 2]) .* precoder, 3);
    % normalize waveform power
    waveform = sqrt(txPower) * waveform / norm(waveform, 'fro');
    % \boldsymbol{X}
    waveformMatrix = vec(waveform) * vec(waveform)';

    % * get block diagonal channel matrix and initialize auxiliary variables
    % \boldsymbol{M}
    [matrixChannel] = matrix_channel(channel);
    % \boldsymbol{t}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        for iSubband = 1 : nSubbands
            auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * waveformMatrix);
        end
    end

    % * update waveform matrix and auxiliary variables iteratively
    isConverged = false;
    counter = 0;
    % \boldsymbol{A}_0
    A0 = - diag(3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        counter = counter + 1;
        % \bar{c}
        cBar = zeros(1, nUsers);
        % \boldsymbol{C}_1
        C1 = cell(1, nUsers);
        % \boldsymbol{A}_1
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
        [~, userIndex] = max(target);

        % * Rank reduction for separable SDP
        waveformMatrix_ = highRankWaveformMatrix;
        waveformRank = rank(waveformMatrix_);
        while waveformRank ^ 2 > nUsers
            % decompose waveform matrix as a product of a matrix V and its Hermitian
            waveformFactor = decompose(waveformMatrix_);

            % % flatten trace equations to standard linear equations
            % coefficient = zeros(nUsers, waveformRank ^ 2);
            % for iUser = 1 : nUsers
            %     coefficient(iUser, :) = reshape((waveformFactor' * (A1{iUser} - A1{userIndex}) * waveformFactor).', 1, []);
            % end
            % % obtain an orthonormal basis for the null space of the coefficient matrix
            % delta = null(coefficient);
            % % nonzero solution can be obtained as a linear combination of null space basis vectors
            % delta = reshape(delta(:, 1), [waveformRank, waveformRank]);
            % % ensure positive semidefiniteness
            % delta = (delta + delta') / 2;

            deltaInit = eye(waveformRank);

            options = optimset('Algorithm', 'Levenberg-Marquardt', 'TolFun', eps, 'TolX', eps, 'Display', 'off', 'MaxIter', 200);
            delta = fsolve(@(delta) rr_equations(delta, waveformFactor, A1, nTxs, nSubbands, nUsers, userIndex), deltaInit, options);

            % calculate eigenvalues of delta
            d = eig(delta);
            % there can be multiple candidates with largest magnitude for each user and we only use the minimum one (the negative one if there is both positive and negative)
            dominantEigenvalue = min(d(abs(d) == max(abs(d))));
            clearvars d;
            waveformMatrix_ = waveformFactor * (eye(waveformRank) - 1 / dominantEigenvalue * delta) * waveformFactor';
            % % ! ensure positive semidefiniteness
            % [v, d] = eig(waveformMatrix_);
            % d(d < 0) = 0;
            % waveformMatrix_ = v * d * v';
            % clearvars v d;
            % update matrix rank
            waveformRank = rank(waveformMatrix_);
        end

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * waveformMatrix_);
            end
        end

        % test convergence
        temp = (norm(waveformMatrix_ - waveformMatrix, 'fro')) / norm(waveformMatrix_, 'fro')
        if (norm(waveformMatrix_ - waveformMatrix, 'fro')) / norm(waveformMatrix_, 'fro') <= tolerance || counter >= 1e2
            isConverged = true;
        end
        waveformMatrix = waveformMatrix_;

    end

    % obtain the rank-1 beamforming vector
    [v, d] = svd(waveformMatrix);
    waveform = v(:, 1) * sqrt(d(1));

    % % decompose waveform matrix for the waveform vector
    % waveform = decompose(waveformMatrix);

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

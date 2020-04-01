function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_che_rr(beta2, beta4, txPower, channel, tolerance, weight, pathloss)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the minimum voltage with rank reduction power allocation based on large-scale fading
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
    %   - for multi-user MISO systems with number of user no larger than 3
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - in each iteration, we first obtain the high rank covariance matrix by CVX, then perform rank reduction for a rank-1 solution
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %   - Y. Huang and D. Palomar, "Rank-Constrained Separable Semidefinite Programming With Applications to Optimal Beamforming," IEEE Transactions on Signal Processing, vol. 58, no. 2, pp. 664–678, 2010.
    %
    % Author & Date: Yang (i@snowztail.com) - 23 Mar 20



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
    A0 = - diag(3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        counter = counter + 1;
        % \bar{c}_q'
        cBar = zeros(nUsers, 1);
        % \boldsymbol{C}'_{q, 1}
        C1 = cell(nUsers, 1);
        % \boldsymbol{A}'_{q, 1}
        A1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            cBar(iUser) = - real(conj(auxiliary(iUser, :)) * A0 * auxiliary(iUser, :).' * nTxs ^ 2 * txPower ^ 2 * pathloss(iUser) ^ 4);
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
                target(iUser) = trace(A1{iUser} * highRankcarrierWeightMatrix(:, :, iUser)) + cBar(iUser);
                traceSum = traceSum + pathloss(iUser) * trace(highRankcarrierWeightMatrix(:, :, iUser));
            end
            minimize(max(target));
            subject to
                traceSum == 1;
        cvx_end
        [~, userIndex] = max(target);

        % * Rank reduction for separable SDP
        carrierWeightMatrix_ = highRankcarrierWeightMatrix;
        carrierWeightRank = zeros(nUsers, 1);
        for iUser = 1 : nUsers
            carrierWeightRank(iUser) = rank(carrierWeightMatrix_(:, :, iUser));
        end
        % carrierWeightRank = sum(carrierWeightRank);
        carrierWeightFactor = cell(nUsers, 1);
        deltaRank = carrierWeightRank;
        while sum(carrierWeightRank .^ 2) > nUsers
            for iUser = 1 : nUsers
                % decompose weight matrix as a product of a matrix V and its Hermitian
                carrierWeightFactor{iUser} = decompose(carrierWeightMatrix_(:, :, iUser));
            end
%             % flatten trace equations to standard linear equations
%             coefficient = zeros(nUsers, sum(carrierWeightRank .^ 2));
%             for iUser = 1 : nUsers
%                 for jUser = 1 : nUsers
%                     coefficient(iUser, 1 + sum(carrierWeightRank(1 : jUser - 1) .^ 2) : sum(carrierWeightRank(1 : jUser) .^ 2)) = reshape((carrierWeightFactor{jUser}' * (A1{iUser} - A1{userIndex}) * carrierWeightFactor{jUser}).', 1, []);
%                 end
%             end
%             % obtain an orthonormal basis for the null space of the coefficient matrix
%             delta = null(coefficient);
%             % nonzero solution can be obtained as a linear combination of null space basis vectors
%             deltaInstance = cell(nUsers, 1);
%             % dominantEigenvalue = zeros(nUsers, 1);
%             deltaRank = zeros(nUsers, 1);
%             for iUser = 1 : nUsers
%                 deltaInstance{iUser} = reshape(delta(1 + sum(carrierWeightRank(1 : iUser - 1) .^ 2) : sum(carrierWeightRank(1 : iUser) .^ 2), 1), [carrierWeightRank(iUser), carrierWeightRank(iUser)]);
%                 % ensure positive semidefiniteness
%                 deltaInstance{iUser} = (deltaInstance{iUser} + deltaInstance{iUser}') / 2;
%                 % prioritize users with high-rank delta
%                 deltaRank(iUser) = rank(deltaInstance{iUser});
%                 % % calculate eigenvalues of delta
%                 % d = eig(deltaInstance{iUser});
%                 % % there can be multiple candidates with largest magnitude for each user and we only use the minimum one (the negative one if there is both positive and negative)
%                 % dominantEigenvalue(iUser) = min(d(abs(d) == max(abs(d))));
%                 % clearvars d;
%             end
%             [~, deltaIndex] = max(deltaRank);
%             % calculate eigenvalues of delta
%             d = eig(deltaInstance{deltaIndex});
%             % there can be multiple candidates with largest magnitude for each user and we only use the minimum one (the negative one if there is both positive and negative)
%             dominantEigenvalue = min(d(abs(d) == max(abs(d))));
%
%             for iUser = 1 : nUsers
%                 % update beamforming matrices
%                 carrierWeightMatrix_(:, :, iUser) = carrierWeightFactor{iUser} * (eye(carrierWeightRank(iUser)) - 1 / dominantEigenvalue * deltaInstance{iUser}) * carrierWeightFactor{iUser}';
%                 % % ! ensure positive semidefiniteness
%                 % [v, d] = eig(carrierWeightMatrix_(:, :, iUser));
%                 % d = real(d);
%                 % d(d < 0) = 0;
%                 % carrierWeightMatrix_(:, :, iUser) = v * d * v';
%                 % carrierWeightMatrix_(:, :, iUser) = (carrierWeightMatrix_(:, :, iUser) + carrierWeightMatrix_(:, :, iUser)') / 2;
%                 % clearvars v d;
%                 carrierWeightRank(iUser) = rank(carrierWeightMatrix_(:, :, iUser));
%             end

%%
            deltaInit = cell(nUsers, 1);
            for iUser = 1 : nUsers
                deltaInit{iUser} = vec(eye(carrierWeightRank(iUser)));
            end
            deltaInit = cat(1, deltaInit{:});

%             deltaInit = cell(nUsers, 1);
%             for iUser = 1 : nUsers
%                 deltaInit{iUser} = (eye(carrierWeightRank(iUser)));
%             end


            % (deltaVector, factor, term, pathloss, nSubbands, nUsers, userIndex, deltaRank)
            % options = optimset('Algorithm', 'levenberg-marquardt', 'Display', 'off');
            options = optimset('Algorithm', 'Levenberg-Marquardt', 'TolFun', eps, 'TolX', eps, 'Display', 'off', 'MaxIter', 200);
            delta = fsolve(@(delta) rr_equations_che(delta, carrierWeightFactor, A1, pathloss, nSubbands, nUsers, userIndex, carrierWeightRank), deltaInit, options);

            % nonzero solution can be obtained as a linear combination of null space basis vectors
            deltaInstance = cell(nUsers, 1);
            % dominantEigenvalue = zeros(nUsers, 1);
            for iUser = 1 : nUsers
                deltaInstance{iUser} = reshape(delta(1 + sum(carrierWeightRank(1 : iUser - 1) .^ 2) : sum(carrierWeightRank(1 : iUser) .^ 2), 1), [carrierWeightRank(iUser), carrierWeightRank(iUser)]);
                % prioritize users with high-rank delta
                deltaRank(iUser) = rank(deltaInstance{iUser});
                % % calculate eigenvalues of delta
                % d = eig(deltaInstance{iUser});
                % % there can be multiple candidates with largest magnitude for each user and we only use the minimum one (the negative one if there is both positive and negative)
                % dominantEigenvalue(iUser) = min(d(abs(d) == max(abs(d))));
                % clearvars d;
            end

            [~, deltaIndex] = max(deltaRank);
            % calculate eigenvalues of delta
            d = eig(deltaInstance{deltaIndex});
            % there can be multiple candidates with largest magnitude for each user and we only use the minimum one (the negative one if there is both positive and negative)
            dominantEigenvalue = min(d(abs(d) == max(abs(d))));

            for iUser = 1 : nUsers
                % update beamforming matrices
                carrierWeightMatrix_(:, :, iUser) = carrierWeightFactor{iUser} * (eye(carrierWeightRank(iUser)) - 1 / dominantEigenvalue * deltaInstance{iUser}) * carrierWeightFactor{iUser}';
                % % ! ensure positive semidefiniteness
                % [v, d] = eig(carrierWeightMatrix_(:, :, iUser));
                % d = real(d);
                % d(d < 0) = 0;
                % carrierWeightMatrix_(:, :, iUser) = v * d * v';
                % carrierWeightMatrix_(:, :, iUser) = (carrierWeightMatrix_(:, :, iUser) + carrierWeightMatrix_(:, :, iUser)') / 2;
                % clearvars v d;
                carrierWeightRank(iUser) = rank(carrierWeightMatrix_(:, :, iUser));
            end
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
            target(iUser) = real(trace(A1{iUser} * carrierWeightMatrix_(:, :, iUser)) + cBar(iUser));
        end
        maxTarget_ = max(target);

        % test convergence
        temp = abs(maxTarget_ - maxTarget) / abs(maxTarget_)
        if abs(maxTarget_ - maxTarget) / abs(maxTarget_) <= tolerance || counter >= 1e2
            isConverged = true;
        end
        carrierWeightMatrix = carrierWeightMatrix_;
        maxTarget = maxTarget_;




        %%




        % Update p_vec_opt_mat

%         for IndUser = 1:NumUsers
%             X_opt_Rank1 = X_opt_RankConstrained_mats(:, :, IndUser);
%             [V_X_opt_Rank1, D_X_opt_Rank1, ~] = svd(X_opt_Rank1);
% %             p_vec_opt_mat(:, IndUser) = V_X_opt_Rank1(:, 1).*sqrt(D_X_opt_Rank1(1,1));
%             p_vec_opt_mat(:, IndUser) = decompose(X_opt_Rank1);
%         end



%         XErrorSqr = norm(carrierWeight_ - carrierWeight, 'fro')^2;
%         XSqr = norm(carrierWeight_, 'fro')^2;
%
%         RelativeError = sqrt(XErrorSqr)/sqrt(XSqr);
%         if RelativeError <= tolerance
%             isConverged = true;
%         end
%
%
%         % test convergence
%         temp = (norm(carrierWeight_ - carrierWeight, 'fro')) / norm(carrierWeight_, 'fro')
%         if (norm(carrierWeight_ - carrierWeight, 'fro')) / norm(carrierWeight_, 'fro') <= tolerance
%             isConverged = true;
%         end
%         carrierWeight = carrierWeight_;

        % test convergence
%         temp = abs(maxTarget_ - maxTarget) / abs(maxTarget_)
%         if abs(maxTarget_ - maxTarget) / abs(maxTarget_) <= tolerance || counter >= 1e2
%             isConverged = true;
%         end
%         carrierWeightMatrix = carrierWeightMatrix_;
%         maxTarget = maxTarget_;

    end
    % obtain the rank-1 beamforming vector
    for iUser = 1 : nUsers
        [v, d] = svd(carrierWeightMatrix(:, :, iUser));
        carrierWeight(:, iUser) = v(:, 1) * sqrt(d(1));
    end

%     % decompose beamforming matrix for the beamforming vector
%     for iUser = 1 : nUsers
%         carrierWeight(:, iUser) = decompose(carrierWeightMatrix(:, :, iUser));
%     end


    % \bar{\boldsymbol{s}}_n
    normalizedWaveform = sum(repmat(reshape(carrierWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* conj(channel), 3) / sqrt(nTxs);
    % \boldsymbol{s}_{\text{asym}}
    waveform = sqrt(txPower) * normalizedWaveform / norm(normalizedWaveform, 'fro');
    % \sum v_{\text{out}}, v\{\text{out}, q}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

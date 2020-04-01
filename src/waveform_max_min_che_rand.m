function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_che_rand(beta2, beta4, txPower, channel, tolerance, pathloss)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the minimum voltage with randomized power allocation based on large-scale fading
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - pathloss [\boldsymbol{\Lambda}] (1 * nUsers): large-scale channel strength reduction
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_{\text{asym}}] (nTxs * nSubbands): the asymptotically optimal complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - maximize the minimum user voltage
    %   - for multi-user MISO systems with arbitrary number of user
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - we first obtain a high rank covariance matrix by CVX, then generate waveform based on one randomized vectors
    %   - only one randomized step is used and the output voltage can vary even for a given channel
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %   - Y. Huang and D. P. Palomar, "Randomized Algorithms for Optimal Solutions of Double-Sided QCQP With Applications in Signal Processing," IEEE Transactions on Signal Processing, vol. 62, no. 5, pp. 1093–1108, 2014.
    %
    % Author & Date: Yang (i@snowztail.com) - 25 Mar 20



    % * initialize complex carrier weight by channel strength
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * permute(sum(channel, 1), [2, 3, 1]) / norm(squeeze(sum(channel)), 'fro');
    % \boldsymbol{X}
    carrierWeightMatrix = zeros(nSubbands, nSubbands, nUsers);
    for iUser = 1 : nUsers
        carrierWeightMatrix(:, :, iUser) = carrierWeight(:, iUser) * carrierWeight(:, iUser)';
    end

    % ! unify notation
    pathloss = 1 ./ pathloss;
    eirp = txPower * nTxs;

    % * get shift matrix and initialize auxiliary variables
    % \boldsymbol{M}'
    [matrixShift] = matrix_shift(channel);
    % \boldsymbol{t}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        for iSubband = 1 : nSubbands
            auxiliary(iUser, iSubband) = trace(matrixShift{iSubband} * carrierWeightMatrix(:, :, iUser));
        end
    end

    % * update carrier weight matrix and auxiliary variables iteratively
    isConverged = false;
    maxTarget = 0;
    % \boldsymbol{A}_0
    A0 = - diag(3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
        % * update term A1', C1', and cBar'
        % \bar{c}'
        cBar = zeros(nUsers, 1);
        % \boldsymbol{C}'_1
        C1 = cell(nUsers, 1);
        % \boldsymbol{A}'_1
        A1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            cBar(iUser) = - real(conj(auxiliary(iUser, :)) * A0 * auxiliary(iUser, :).' * eirp ^ 2 * pathloss(iUser) ^ 4);
            C1{iUser} = - (beta2 * eirp * pathloss(iUser) ^ 2 + 3 * eirp ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * matrixShift{1};
            if nSubbands > 1
                C1{iUser} = C1{iUser} - 3 * beta4 * eirp ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, matrixShift{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end

        % * solve high rank carrier weight matrix in SDP problem by cvx
        cvx_begin quiet
            % \boldsymbol{X}
            variable solutionMatrix(nSubbands, nSubbands, nUsers) complex semidefinite;
            % Tr{\boldsymbol{AX}} + \bar{C}
            target = cvx(zeros(1, nUsers));
            traceSum = 0;
            for iUser = 1 : nUsers
                target(iUser) = trace(A1{iUser} * solutionMatrix(:, :, iUser)) + cBar(iUser);
                traceSum = traceSum + pathloss(iUser) * trace(solutionMatrix(:, :, iUser));
            end
            minimize(max(target));
            subject to
                traceSum == 1;
        cvx_end
        [~, userIndex] = max(target);

        % * update term B1, B2, and Q
        B1 = A1;
        B1{userIndex} = 0;
        B2 = cell(nUsers, 1);
        Q = cell(nUsers, 1);
        carrierWeight_ = zeros(nSubbands, nUsers);
        carrierWeightMatrix_ = zeros(nSubbands, nSubbands, nUsers);
        for iUser = 1 : nUsers
            B2{iUser} = pathloss(iUser) * eye(nSubbands);
            [U, ~] = eig(solutionMatrix(:, :, iUser) ^ (1 / 2) * B1{iUser} * solutionMatrix(:, :, iUser) ^ (1 / 2));
            Q{iUser} = U' * solutionMatrix(:, :, iUser) ^ (1 / 2) * B2{iUser} * solutionMatrix(:, :, iUser) ^ (1 / 2) * U;
            % * input generated matrix into the randomization algorithm
            randomizedVector = randomized_solution(Q{iUser});
            carrierWeight_(:, iUser) = solutionMatrix(:, :, iUser) ^ (1 / 2) * U * randomizedVector;
            carrierWeightMatrix_(:, :, iUser) = carrierWeight(:, iUser) * carrierWeight(:, iUser)';
            clearvars U;
        end
        % Update - \gamma_2 and \boldsymbol{t}
        for iUser = 1 : nUsers
            target(iUser) = real(trace(A1{iUser} * carrierWeightMatrix_(:, :, iUser)) + cBar(iUser));
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixShift{iSubband} * carrierWeightMatrix_(:, :, iUser));
            end
        end
        maxTarget_ = max(target);

        % * test convergence
        if abs(maxTarget_ - maxTarget) / abs(maxTarget_) <= tolerance
            isConverged = true;
        end
        carrierWeight = carrierWeight_;
        maxTarget = maxTarget_;
    end

    % * the asymptotic optimal spatial precoder is asymptotic matched filter (divide by nTxs rather than channel norm)
    % \bar{\tilde{s}}
    precoder = conj(channel) / sqrt(nTxs);

    % * construct waveform
    % \boldsymbol{s}
    waveform = sum(permute(carrierWeight, [3, 1, 2]) .* precoder, 3);
    % normalize waveform power
    waveform = sqrt(txPower) * waveform / norm(waveform, 'fro');

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

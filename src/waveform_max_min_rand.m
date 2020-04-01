function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_max_min_rand(beta2, beta4, txPower, channel, tolerance, weight, nCandidates)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the minimum voltage with randomized power allocation
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w] (1 * nUsers): user weights
    %   - nCandidates: number of random feasible rank-1 solutions to generate
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
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
    % \boldsymbol{A}_0
    A0 = - diag(3 * beta4 * [1 / 2, ones(1, nSubbands - 1)]);
    while ~isConverged
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
                C1{iUser} = C1{iUser} - 3 * beta4 * sum(cat(3, matrixChannel{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end

        % * solve high rank waveform matrix in SDP problem by cvx
        cvx_begin quiet
            % \boldsymbol{X}
            variable solutionMatrix(nTxs * nSubbands, nTxs * nSubbands) complex semidefinite;
            % Tr{\boldsymbol{AX}} + \bar{C}
            target = cvx(zeros(1, nUsers));
            for iUser = 1 : nUsers
                target(iUser) = trace(A1{iUser} * solutionMatrix) + cBar(iUser);
            end
            minimize(max(target));
            subject to
                trace(solutionMatrix) <= txPower;
        cvx_end
        waveformMatrix_ = solutionMatrix;
        % Update \boldsymbol{t}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * waveformMatrix_);
            end
        end

        % * test convergence
        if (norm(waveformMatrix_ - waveformMatrix, 'fro')) / norm(waveformMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        waveformMatrix = waveformMatrix_;
    end

    % * choose an optimal waveform constructed from random vectors
    % \hat{\boldsymbol{x}}
    candidateWaveform = cell(nCandidates, 1);
    % - \gamma_2
    candidateTarget = zeros(nCandidates, nUsers);
    [v, d] = eig(waveformMatrix);
    for iCandidate = 1 : nCandidates
        candidateWaveform{iCandidate} = v * d ^ (1 / 2) * exp(1i * 2 * pi * rand(nTxs * nSubbands, 1));
        % \boldsymbol{X}
        candidateWaveformMatrix = candidateWaveform{iCandidate} * candidateWaveform{iCandidate}';
        % Update \boldsymbol{t}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannel{iUser, iSubband} * candidateWaveformMatrix);
            end
        end
        % update related terms
        for iUser = 1 : nUsers
            cBar(iUser) = - real(conj(auxiliary(iUser, :)) * A0 * auxiliary(iUser, :).');
            C1{iUser} = - (beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * matrixChannel{iUser, 1};
            if nSubbands > 1
                C1{iUser} = C1{iUser} - weight(iUser) * 3 * beta4 * sum(cat(3, matrixChannel{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = C1{iUser} + C1{iUser}';
        end
        for iUser = 1 : nUsers
            candidateTarget(iCandidate, iUser) = real(trace(A1{iUser} * candidateWaveformMatrix) + cBar(iUser));
        end
    end
    [~, candidateIndex] = min(max(candidateTarget, [], 2));
    waveform = reshape(candidateWaveform{candidateIndex}, [nTxs, nSubbands]);

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

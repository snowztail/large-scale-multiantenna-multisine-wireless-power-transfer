function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_wsums(beta2, beta4, txPower, channel, tolerance, weight)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the weighted sum voltage by using suboptimal precoder and optimizing carrier weight
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [\boldsymbol{w}] (1 * nUsers): user weights
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
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
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20



    % * initialize complex carrier weight (thus power allocation) by channel strength
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * permute(sum(channel, 1), [2, 3, 1]) / norm(squeeze(sum(channel)), 'fro');
    % \boldsymbol{X}
    carrierWeightMatrix = carrierWeight * carrierWeight';

    % * considering the linear model, the optimal spatial beamformer is the dominant eigenvector of the weighted channel matrix
    % \boldsymbol{w}
    precoder = zeros(nTxs, nSubbands);
    for iSubband = 1 : nSubbands
        [v, d] = eig(conj(squeeze(channel(:, iSubband, :))) * diag(weight) * squeeze(channel(:, iSubband, :)).');
        precoder(:, iSubband) = v(:, diag(d) == max(diag(d)));
        clearvars v d;
    end
    % unify notations
    precoder = repmat(precoder, [1, 1, nUsers]);

    % * update equivalent channel matrix and auxiliary variables iteratively
    % \boldsymbol{M}'''
    [matrixChannelEquivalent] = matrix_channel_equivalent(channel, precoder);
    % \boldsymbol{t}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        for iSubband = 1 : nSubbands
            auxiliary(iUser, iSubband) = trace(matrixChannelEquivalent{iUser, iSubband} * carrierWeightMatrix);
        end
    end

    % * update carrier weight matrix and auxiliary variables iteratively
    isConverged = false;
    while ~isConverged
        % * update term A1''' and C1'''
        % \boldsymbol{C}'''_1
        C1 = 0;
        for iUser = 1 : nUsers
            C1 = C1 - weight(iUser) * ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * matrixChannelEquivalent{iUser, 1});
            if nSubbands > 1
                C1 = C1 - weight(iUser) * (3 * beta4 * sum(cat(3, matrixChannelEquivalent{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3));
            end
        end
        % \boldsymbol{A}'''_1
        A1 = C1 + C1';

        % * solve rank-1 carrier weight matrix in closed form
        [v, d] = eig(A1);
        % \boldsymbol{x}^{\star}
        carrierWeight_ = sqrt(txPower) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        carrierWeightMatrix_ = carrierWeight_ * carrierWeight_';
        % Update \boldsymbol{t}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(matrixChannelEquivalent{iUser, iSubband} * carrierWeightMatrix_);
            end
        end

        % * test convergence
        if (norm(carrierWeightMatrix_ - carrierWeightMatrix, 'fro')) / norm(carrierWeightMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        carrierWeight = carrierWeight_;
        carrierWeightMatrix = carrierWeightMatrix_;
    end

    % * construct waveform
    % \boldsymbol{s}
    waveform = sum(permute(carrierWeight, [3, 1, 2]) .* precoder, 3);
    % normalize waveform power
    waveform = sqrt(txPower) * waveform / norm(waveform, 'fro');

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

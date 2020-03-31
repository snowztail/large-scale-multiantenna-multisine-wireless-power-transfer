function [waveform, sumVoltage, userVoltage] = waveform_wsums(beta2, beta4, txPower, channel, tolerance, weight)
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
    % Author & Date: Yang (i@snowztail.com) - 11 Mar 20



    [nTxs, nSubbands, nUsers] = size(channel);
    % ? initialize \boldsymbol{p} by uniform power allocation
    carrierWeight = sqrt(txPower / nSubbands) * ones(nSubbands, 1);
    % \boldsymbol{X}
    carrierWeightMatrix = carrierWeight * carrierWeight';

    % \boldsymbol{w}_n
    precoder = zeros(nTxs, nSubbands);
    for iSubband = 1 : nSubbands
        % the optimal spatial beamformer is the dominant eigenvector of the weighted channel matrix
        [v, d] = eig(conj(squeeze(channel(:, iSubband, :))) * diag(weight) * squeeze(channel(:, iSubband, :)).');
        precoder(:, iSubband) = v(:, diag(d) == max(diag(d)));
        clearvars v d;
    end

    % \boldsymbol{M}'''_{q, k}
    equivalentChannelMatrix = cell(nUsers, nSubbands);
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        % \boldsymbol{h}_{e, q}
        equivalentSubchannel = diag(precoder' * conj(channel(:, :, iUser)));
        % \boldsymbol{M}'''_q
        equivalentSubchannelMatrix = equivalentSubchannel * equivalentSubchannel';
        for iSubband = 1 : nSubbands
            equivalentChannelMatrix{iUser, iSubband} = diag(diag(equivalentSubchannelMatrix, iSubband - 1), iSubband - 1);
            auxiliary(iUser, iSubband) = trace(equivalentChannelMatrix{iUser, iSubband} * carrierWeightMatrix);
        end
    end

    isConverged = false;
    while ~isConverged
        % \boldsymbol{C}'''_1
        C1 = 0;
        for iUser = 1 : nUsers
            C1 = C1 - weight(iUser) * ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * equivalentChannelMatrix{iUser, 1});
            if nSubbands > 1
                C1 = C1 - weight(iUser) * (3 * beta4 * sum(cat(3, equivalentChannelMatrix{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3));
            end
        end
        % \boldsymbol{A}'''_1
        A1 = C1 + C1';

        % * Solve rank-1 \boldsymbol{X}^{\star} in closed form (low complexity)
        % \boldsymbol{x}^{\star}
        [v, d] = eig(A1);
        carrierWeight = sqrt(txPower) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        carrierWeightMatrix_ = carrierWeight * carrierWeight';
        % Update \boldsymbol{t}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(equivalentChannelMatrix{iUser, iSubband} * carrierWeightMatrix_);
            end
        end
        % test convergence
        if (norm(carrierWeightMatrix_ - carrierWeightMatrix, 'fro')) / norm(carrierWeightMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        carrierWeightMatrix = carrierWeightMatrix_;
    end

    % v_{\text{out}, q}
    userVoltage = zeros(1, nUsers);
    for iUser = 1 : nUsers
        userVoltage(iUser) = beta2 * carrierWeight' * equivalentChannelMatrix{iUser, 1} * carrierWeight + (3 / 2) * beta4 * carrierWeight' * equivalentChannelMatrix{iUser, 1} * carrierWeight * (carrierWeight' * equivalentChannelMatrix{iUser, 1} * carrierWeight)';
        if nSubbands > 1
            for iSubband = 1 : nSubbands - 1
                userVoltage(iUser) = userVoltage(iUser) + 3 * beta4 * carrierWeight' * equivalentChannelMatrix{iUser, iSubband + 1} * carrierWeight * (carrierWeight' * equivalentChannelMatrix{iUser, iSubband + 1} * carrierWeight)';
            end
        end
    end
    userVoltage = real(userVoltage);
    % \boldsymbol{s_n}
    waveform = carrierWeight.' .* precoder;
    % \sum v_{\text{out}}
    sumVoltage = sum(userVoltage);

end

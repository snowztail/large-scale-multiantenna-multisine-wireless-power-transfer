function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_che_wsum(beta2, beta4, txPower, channel, tolerance, weight, pathloss)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize weighted sum voltage with power allocation based on large-scale fading
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w] (1 * nUsers): user weights
    %   - pathloss [\Lambda] (1 * nUsers): user pathlosses
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_{\text{asym}}] (nTxs * nSubbands): the asymptotically optimal complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - for multi-user MISO systems
    %   - exploit channel hardening
    %   - consider path loss in power allocation and assume fading is i.i.d. in space and frequency
    %   - the spatial beamformer is still a function of short-term CSI
    %   - note the notation difference with the single-user case
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 12 Mar 20



    % * initialize complex carrier weight by channel strength
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * permute(sum(channel, 1), [2, 3, 1]) / norm(squeeze(sum(channel)), 'fro');
    % \boldsymbol{X}
    carrierWeightMatrix = carrierWeight * carrierWeight';

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
            auxiliary(iUser, iSubband) = carrierWeight(:, iUser)' * matrixShift{iSubband} * carrierWeight(:, iUser);
        end
    end

    % * update carrier weight matrix and auxiliary variables iteratively
    isConverged = false;
    while ~isConverged
        % * update term A1' and C1'
        % \boldsymbol{C}'_1
        C1 = cell(nUsers, 1);
        % \boldsymbol{A}'_1
        A1 = cell(nUsers, 1);
        for iUser = 1 : nUsers
            C1{iUser} = - ((beta2 * eirp * pathloss(iUser) ^ 2 + 3 * eirp ^ 2 * pathloss(iUser) ^ 4 * beta4 * auxiliary(iUser, 1)) / 2 * matrixShift{1});
            if nSubbands > 1
                C1{iUser} = C1{iUser} - 3 * beta4 * eirp ^ 2 * pathloss(iUser) ^ 4 * sum(cat(3, matrixShift{2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
            A1{iUser} = weight(iUser) * (C1{iUser} + C1{iUser}');
        end
        A1 = blkdiag(A1{:});
        % \boldsymbol{\Lambda}
        pathlossMatrix = diag(repelem(pathloss, nSubbands));

        % * solve rank-1 carrier weight matrix in closed form
        [v, d] = eig(pathlossMatrix \ A1);
        vMin = v(:, diag(d) == min(diag(d)));
        % \bar{\boldsymbol{p}}
        carrierWeight_ = sqrt(1 / (vMin' * pathlossMatrix * vMin)) * vMin;
        clearvars v d vMin;
        % \boldsymbol{p}^{\star}
        carrierWeight_ = reshape(carrierWeight_, nSubbands, nUsers);
        % \boldsymbol{X}^{\star}
        carrierWeightMatrix_ = carrierWeight_ * carrierWeight_';
        % update \boldsymbol{t}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = carrierWeight_(:, iUser)' * matrixShift{iSubband} * carrierWeight_(:, iUser);
            end
        end

        % * test convergence
        if (norm(carrierWeightMatrix_ - carrierWeightMatrix, 'fro')) / norm(carrierWeightMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        carrierWeight = carrierWeight_;
        carrierWeightMatrix = carrierWeightMatrix_;
    end

    % * the asymptotic optimal spatial precoder is asymptotic matched filter (divide by nTxs rather than channel norm)
    % \bar{\tilde{s}}
    precoder = conj(channel) / sqrt(nTxs);

    % * construct waveform
    % \bar{\boldsymbol{s}}
    waveform = sum(permute(carrierWeight, [3, 1, 2]) .* precoder, 3);
    % \boldsymbol{s}_{\text{asym}}
    waveform = sqrt(txPower) * waveform / norm(waveform, 'fro');

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

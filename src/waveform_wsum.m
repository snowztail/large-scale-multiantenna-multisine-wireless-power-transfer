function [waveform, sumVoltage, userVoltage, minVoltage] = waveform_wsum(beta2, beta4, txPower, channel, tolerance, weight)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %   - maximize the weight sum volgate by designing waveform directly
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
    %   - waveform [\boldsymbol{s}] (nTxs * nSubbands * nUsers): complex waveform weights for each transmit antenna and subband
    %   - sumVoltage [\sum v_{\text{out}}]: sum of rectifier output DC voltage over all users
    %   - userVoltage [v_{\text{out}}]: individual user voltages
    %   - minVoltage [\min v_{\text{out}}]: minimum user voltage
    %
    % Comment(s):
    %   - for multi-user MISO systems
    %   - the frequency domain power allocation and spatial beamforming are coupled for multi-user scenario
    %   - optimize frequency domain power allocation and spatial beamforming jointly
    %   - note the notation difference with the single-user case
    %   - obtain the rank-1 frequency weight matrix in closed form for a complexity of \mathcal{O}((MN) ^ 3)
    %   - as all eigenvalues of \boldsymbol{A}_1 are negative, we choose the eigenvector corresponding to the minimum eigenvalue (maximum in magnitude)
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 10 Mar 20



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
    while ~isConverged
        % * update term A1 and C1
        % \boldsymbol{C}_1
        C1 = 0;
        for iUser = 1 : nUsers
            C1 = C1 - weight(iUser) * ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * matrixChannel{iUser, 1});
            if nSubbands > 1
                C1 = C1 - weight(iUser) * 3 * beta4 * sum(cat(3, matrixChannel{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser, 2 : end)), [1, 1, nSubbands - 1]), 3);
            end
        end
        % \boldsymbol{A}_1
        A1 = C1 + C1';

        % * solve rank-1 waveform in closed form
        [v, d] = eig(A1);
        % \boldsymbol{x}^{\star}
        waveform_ = sqrt(txPower) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        waveformMatrix_ = waveform_ * waveform_';
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
        waveform = waveform_;
        waveformMatrix = waveformMatrix_;
    end

    % * construct waveform
    % \boldsymbol{s}
    waveform = reshape(waveform, [nTxs, nSubbands]);

    % * compute output voltages
    % \sum v_{\text{out}}, v\{\text{out}}, \min v_{\text{out}}
    [sumVoltage, userVoltage, minVoltage] = harvester_compact(beta2, beta4, waveform, channel);

end

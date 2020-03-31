function [waveform, voltage] = waveform_su(beta2, beta4, txPower, channel, tolerance)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - txPower [P]: transmit power constraint
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s}_n] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - voltage [\sum v_{\text{out}}]: rectifier output DC voltage
    %
    % Comment(s):
    %   - for single-user MISO systems
    %   - the optimal spatial beamforming for SU WPT is MRT
    %   - obtain the rank-1 frequency weight matrix in closed form for a complexity of \mathcal{O}(N ^ 3)
    %   - as all eigenvalues of \boldsymbol{A}''_1 are negative, we choose the eigenvector corresponding to the minimum eigenvalue (maximum in magnitude)
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 06 Mar 20



    % * initialize complex carrier weight (thus power allocation) by channel strength
    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{p}
    carrierWeight = sqrt(txPower) * permute(sum(channel, 1), [2, 3, 1]) / norm(squeeze(sum(channel)), 'fro');
    % \boldsymbol{X}
    carrierWeightMatrix = carrierWeight * carrierWeight';

    % * get channel norm matrix and initialize auxiliary variables
    % \boldsymbol{M}''
    [matrixChannelNorm] = matrix_channel_norm(channel);
    % \boldsynbol{t}
    auxiliary = zeros(1, nSubbands);
    for iSubband = 1 : nSubbands
        auxiliary(iSubband) = trace(matrixChannelNorm{iSubband} * carrierWeightMatrix);
    end

    % * update carrier weight matrix and auxiliary variables iteratively
    isConverged = false;
    while ~isConverged
        % * update term A1'' and C1''
        % \boldsymbol{C}''_1
        C1 = - ((beta2 + 3 * beta4 * auxiliary(1)) / 2 * matrixChannelNorm{1});
        if nSubbands > 1
            C1 = C1 - (3 * beta4 * sum(cat(3, matrixChannelNorm{2 : end}) .* reshape(conj(auxiliary(2 : end)), [1, 1, nSubbands - 1]), 3));
        end
        % \boldsymbol{A}''_1
        A1 = C1 + C1';

        % * solve rank-1 carrier weight matrix in closed form
        [v, d] = eig(A1);
        % \boldsymbol{x}^{\star}
        carrierWeight_ = sqrt(txPower) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        carrierWeightMatrix_ = carrierWeight_ * carrierWeight_';
        % Update \boldsymbol{t}
        for iSubband = 1 : nSubbands
            auxiliary(iSubband) = trace(matrixChannelNorm{iSubband} * carrierWeightMatrix_);
        end

        % * test convergence
        if (norm(carrierWeightMatrix_ - carrierWeightMatrix, 'fro')) / norm(carrierWeightMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        carrierWeight = carrierWeight_;
        carrierWeightMatrix = carrierWeightMatrix_;
    end

    % * optimum single-user precoder is MRT
    % \boldsymbol{\tilde{s}}
    precoder = conj(channel) ./ vecnorm(channel, 2, 1);

    % * construct waveform
    % \boldsymbol{s}
    waveform = sum(repmat(reshape(carrierWeight, [1 nSubbands nUsers]), [nTxs 1 1]) .* precoder, 3);

    % * compute output voltage
    % v_{\text{out}}
    [voltage] = harvester_compact(beta2, beta4, waveform, channel);

end

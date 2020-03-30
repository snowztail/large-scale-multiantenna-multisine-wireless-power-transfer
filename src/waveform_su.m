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



    % * initialize complex carrier weight (thus power allocation) by matched filter
    [~, nSubbands] = size(channel);
    % \xi_n
    % carrierWeight = sqrt(txPower) * conj(channel) / norm(channel, 'fro');
    carrierWeight = sqrt(txPower / nSubbands) * ones(nSubbands, 1);
    % \boldsymbol{X}
    carrierWeightMatrix = carrierWeight * carrierWeight';
    % \boldsymbol{M}''_k
    channelNormMatrix = cell(1, nSubbands);
    % t_k
    auxiliary = zeros(1, nSubbands);
    for iSubband = 1 : nSubbands
        channelNormMatrix{iSubband} = diag(diag(vecnorm(channel, 2, 1).' * vecnorm(channel, 2, 1), iSubband - 1), iSubband - 1);
        auxiliary(iSubband) = trace(channelNormMatrix{iSubband} * carrierWeightMatrix);
    end

    isConverged = false;
    while ~isConverged
        % \boldsymbol{C}''_1
        termC1 = - ((beta2 + 3 * beta4 * auxiliary(1)) / 2 * channelNormMatrix{1});
        if nSubbands > 1
            termC1 = termC1 - (3 * beta4 * sum(cat(3, channelNormMatrix{2 : end}) .* reshape(conj(auxiliary(2 : end)), [1, 1, nSubbands - 1]), 3));
        end
        % \boldsymbol{A}''_1
        termA1 = termC1 + termC1';

        % * Solve rank-1 \boldsymbol{X}^{\star} in closed form (low complexity)
        % \boldsymbol{x}^{\star}
        [v, d] = eig(termA1);
        carrierWeight = sqrt(txPower) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        carrierWeightMatrix_ = carrierWeight * carrierWeight';

        % % * Solve high rank \boldsymbol{X} in SDP problem by cvx (high complexity)
        % cvx_begin
        %     cvx_solver mosek
        %     variable carrierWeightMatrix_
        %     minimize trace(termA1 * carrierWeightMatrix_)
        %     subject to
        %         trace(carrierWeightMatrix_) <= txPower;
        % cvx_end

        % Update \boldsymbol{t}
        for iSubband = 1 : nSubbands
            auxiliary(iSubband) = trace(channelNormMatrix{iSubband} * carrierWeightMatrix_);
        end
        % test convergence
        if (norm(carrierWeightMatrix_ - carrierWeightMatrix, 'fro')) / norm(carrierWeightMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        carrierWeightMatrix = carrierWeightMatrix_;
    end
    % \boldsymbol{\tilde{s}_n}
    spatialPrecoder = conj(channel) ./ vecnorm(channel, 2, 1);
    % \boldsymbol{s_n}
    waveform = carrierWeight.' .* spatialPrecoder;
    % v_{\text{out},q}
    voltage = beta2 * carrierWeight' * channelNormMatrix{1} * carrierWeight + (3 / 2) * beta4 * carrierWeight' * channelNormMatrix{1} * carrierWeight * (carrierWeight' * channelNormMatrix{1} * carrierWeight)';
    if nSubbands > 1
        for iSubband = 1 : nSubbands - 1
            voltage = voltage + 3 * beta4 * carrierWeight' * channelNormMatrix{iSubband + 1} * carrierWeight * (carrierWeight' * channelNormMatrix{iSubband + 1} * carrierWeight)';
        end
    end

end

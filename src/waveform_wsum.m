function [waveform] = waveform_wsum(beta2, beta4, powerBudget, channel, tolerance, weight)
    % Function:
    %   - optimize the amplitude and phase of transmit multisine waveform
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - powerBudget [P]: transmit power constraint
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %   - tolerance [\epsilon]: convergence ratio
    %   - weight [w_q]: user weights
    %
    % OutputArg(s):
    %   - waveform [\boldsymbol{s_n}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
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


    % single receive antenna
    [nTxs, nSubbands, nUsers] = size(channel);
    % ? initialize \boldsymbol{s} by uniform power allocation and omidirectional beamforming
    waveform = sqrt(powerBudget / nTxs / nSubbands) * ones(nTxs, nSubbands);
    % \boldsymbol{X}
    waveformMatrix = waveform(:) * waveform(:)';

    % \boldsymbol{M}_{q, k}
    channelMatrix = cell(nUsers, nSubbands);
    % t_{q, k}
    auxiliary = zeros(nUsers, nSubbands);
    for iUser = 1 : nUsers
        subchannel = channel(:, :, iUser);
        % \boldsymbol{M}_{q}
        subchannelMatrix = conj(subchannel(:)) * subchannel(:).';
        for iSubband = 1 : nSubbands
            channelMatrix{iUser, iSubband} = zeros(nTxs * nSubbands);
            for jSubband = 1 : nSubbands + 1 - iSubband
                channelMatrix{iUser, iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = subchannelMatrix((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            end
            auxiliary(iUser, iSubband) = trace(channelMatrix{iUser, iSubband} * waveformMatrix);
        end
    end

    isConverged = false;
    while ~isConverged
        % \boldsymbol{C}_1
        termC1 = 0;
        for iUser = 1 : nUsers
            termC1 = termC1 - weight(iUser) * ((beta2 + 3 * beta4 * auxiliary(iUser, 1)) / 2 * channelMatrix{iUser, 1} + 3 * beta4 * sum(cat(3, channelMatrix{iUser, 2 : end}) .* reshape(conj(auxiliary(iUser,2 : end)), [1, 1, nSubbands - 1]), 3));
        end
        % \boldsymbol{A}_1
        termA1 = termC1 + termC1';

        % * Solve rank-1 \boldsymbol{X}^{\star} in closed form (low complexity)
        % \boldsymbol{x}^{\star}
        [v, d] = eig(termA1);
        waveform = sqrt(powerBudget) * v(:, diag(d) == min(diag(d)));
        clearvars v d;
        % \boldsymbol{X}^{\star}
        waveformMatrix_ = waveform * waveform';

        % Update \boldsymbol{t}_{q, k}
        for iUser = 1 : nUsers
            for iSubband = 1 : nSubbands
                auxiliary(iUser, iSubband) = trace(channelMatrix{iUser, iSubband} * waveformMatrix_);
            end
        end

        % test convergence
        if (norm(waveformMatrix_ - waveformMatrix, 'fro')) / norm(waveformMatrix_, 'fro') <= tolerance
            isConverged = true;
        end
        waveformMatrix = waveformMatrix_;
    end
    % \boldsymbol{s_n}
    waveform = reshape(waveform, [nTxs, nSubbands]);

end

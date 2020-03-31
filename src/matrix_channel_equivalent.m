function [matrixChannelEquivalent] = matrix_channel_equivalent(channel, precoder)
    % Function:
    %   - calculate the covariance matrix of equivalent channel
    %
    % InputArg(s):
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %   - precoder [\boldsymbol{w}] (nTxs * nSubbands * nUsers): spatial beamformer
    %
    % OutputArg(s):
    %   - matrixChannelEquivalent [\boldsymbol{M}'''] (nUsers * nSubbands) {nSubbands * nSubbands}: covariance matrix of equivalent channel
    %
    % Comment(s):
    %   - the index is above the main diagonal
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Mar 20



    [~, nSubbands, nUsers] = size(channel);
    % \boldsymbol{M}'''
    matrixChannelEquivalent = cell(nUsers, nSubbands);
    for iUser = 1 : nUsers
        % \boldsymbol{h}_e
        subchannel = diag(precoder(:, :, iUser)' * conj(channel(:, :, iUser)));
        for iSubband = 1 : nSubbands
            matrixChannelEquivalent{iUser, iSubband} = diag(diag(subchannel * subchannel', iSubband - 1), iSubband - 1);
        end
    end

end

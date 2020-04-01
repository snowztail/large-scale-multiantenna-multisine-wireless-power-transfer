function [matrixChannel] = matrix_channel(channel)
    % Function:
    %   - calculate the covariance matrix of channel
    %
    % InputArg(s):
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - matrixChannel [\boldsymbol{M}] (nUsers * nSubbands) {(nTxs * nSubbands) * (nTxs * nSubbands)}: only keep the iSubband-th block diagonal in the channel covariance matrix and set the others to 0
    %
    % Comment(s):
    %   - the index is above the main diagonal
    %   - each block is of size nTxs * nTxs
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Mar 20


    [nTxs, nSubbands, nUsers] = size(channel);
    % \boldsymbol{M}
    matrixChannel = cell(nUsers, nSubbands);
    for iUser = 1 : nUsers
        matrixSubchannel = conj(vec(channel(:, :, iUser))) * vec(channel(:, :, iUser)).';
        for iSubband = 1 : nSubbands
            matrixChannel{iUser, iSubband} = zeros(nTxs * nSubbands);
            for jSubband = 1 : nSubbands + 1 - iSubband
                matrixChannel{iUser, iSubband}((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs) = matrixSubchannel((jSubband - 1) * nTxs + 1 : jSubband * nTxs, (iSubband - 1) * nTxs + (jSubband - 1) * nTxs + 1 : (iSubband - 1) * nTxs + jSubband * nTxs);
            end
        end
    end

end

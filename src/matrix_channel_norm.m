function [matrixChannelNorm] = matrix_channel_norm(channel)
    % Function:
    %   - calculate the covariance matrix of channel norm
    %
    % InputArg(s):
    %   - channel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %   - nSubbands [N]: number of subbands/subcarriers
    %
    % OutputArg(s):
    %   - matrixChannelNorm [\boldsymbol{M}''] (nSubbands * 1) {nSubbands * nSubbands}: covariance matrix of channel norm
    %
    % Comment(s):
    %   - for single-user MISO systems
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 30 Mar 20



    [~, nSubbands, ~] = size(channel);
    % \boldsymbol{M}''
    matrixChannelNorm = cell(nSubbands, 1);
    for iSubband = 1 : nSubbands
        matrixChannelNorm{iSubband} = diag(diag(vecnorm(channel, 2, 1).' * vecnorm(channel, 2, 1), iSubband - 1), iSubband - 1);
    end

end

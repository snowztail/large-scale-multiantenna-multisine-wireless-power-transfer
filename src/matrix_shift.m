function [matrixShift] = matrix_shift(channel)
    % Function:
    %   - calculate the shift matrix
    %
    % InputArg(s):
    %   - channel [\boldsymbol{h}] (nTxs * nSubbands * nUsers): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - matrixShift [\boldsymbol{M}'] (1 * nSubbands) { nSubbands * nSubbands}: the iSubband-th diagonal is all ones while other entries are all zeros
    %
    % Comment(s):
    %   - the index is above the main diagonal
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 31 Mar 20



    [~, nSubbands, ~] = size(channel);
    % \boldsymbol{M}'
    matrixShift = cell(1, nSubbands);
    for iSubband = 1 : nSubbands
        matrixShift{iSubband} = diag(diag(ones(nSubbands), iSubband - 1), iSubband - 1);
    end

end

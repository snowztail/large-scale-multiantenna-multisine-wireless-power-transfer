function voltage = harvester(beta2, beta4, waveform, subchannel)
    % Function:
    %   - calculate the harvester output voltage
    %
    % InputArg(s):
    %   - beta2 [\beta_2]: diode second-order parameter
    %   - beta4 [\beta_4]: diode fourth-order parameter
    %   - waveform [\boldsymbol{s_n}] (nTxs * nSubbands): complex waveform weights for each transmit antenna and subband
    %   - subchannel [\boldsymbol{h_{q, n}}] (nTxs * nSubbands): channel frequency response at each subband
    %
    % OutputArg(s):
    %   - voltage [v_{\text{out}}]: rectifier output DC voltage
    %
    % Comment(s):
    %   - truncate the voltage expression to the fourth order to capture fundamental behavior of rectifier nonlinearity
    %
    % Reference(s):
    %   - Y. Huang and B. Clerckx, "Large-Scale Multiantenna Multisine Wireless Power Transfer," IEEE Transactions on Signal Processing, vol. 65, no. 21, pp. 5812–5827, Jan. 2017.
    %
    % Author & Date: Yang (i@snowztail.com) - 08 Mar 20


    % single receive antenna
    [~, nSubbands] = size(subchannel);

    % the second order term in Taylor expansion
    term2 = 0;
    for iSubband = 1 : nSubbands
        term2 = term2 + waveform(:, iSubband)' * conj(subchannel(:, iSubband)) * subchannel(:, iSubband).' * waveform(:, iSubband);
    end
    % the fourth order term in Taylor expansion
    term4 = 0;
    for iSubband1 = 1 : nSubbands
        for iSubband2 = 1 : nSubbands
            for iSubband3 = 1 : nSubbands
                for iSubband4 = 1 : nSubbands
                    % output DC voltage if balanced
                    isBalanced = iSubband1 + iSubband2 == iSubband3 + iSubband4;
                    if isBalanced
                        term4 = term4 + (3 / 2) * waveform(:, iSubband3)' * conj(subchannel(:, iSubband3)) * subchannel(:, iSubband1).' * waveform(:, iSubband1) * waveform(:, iSubband4)' * conj(subchannel(:, iSubband4)) * subchannel(:, iSubband2).' * waveform(:, iSubband2);
                    end
                end
            end
        end
    end
    voltage = real(beta2 * term2 + beta4 * term4);

end

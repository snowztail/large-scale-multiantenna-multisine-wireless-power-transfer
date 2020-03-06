powerBudget = 1;
nTxs = 4;
nSubbands = 8;
beta2 = 2e-3;
beta4 = 4e-3;
tolerance = 1e-3;

% \boldsymbol{h}_q
channel = sqrt(1 / 2) * (randn(nTxs * nSubbands, 1) + 1i * randn(nTxs * nSubbands, 1));

% \boldsymbol{h}_{q,n}
subchannel = reshape(channel, [nTxs, nSubbands]);

frequencyPrecoder = wpt_su(beta2, beta4, powerBudget, nSubbands, nTxs, subchannel, tolerance);

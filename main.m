powerBudget = 1;
nTxs = 4;
nSubbands = 8;
beta4 = 1e-3;

% \boldsymbol{h}_q
channel = randn(nTxs * nSubbands, 1);

% \boldsymbol{h}_{q,n}
subchannel = reshape(channel, [nTxs, nSubbands]);

frequencyPrecoder = wpt_su(beta4, powerBudget, nSubbands, nTxs, subchannel);

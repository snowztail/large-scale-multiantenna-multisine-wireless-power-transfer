%% Transceiver
eirp = db2pow(36 - 30);
nTxs = 4;
powerBudget = eirp / nTxs;
%% Channel
centerFrequency = 2.4e9;
bandwidth = 1e7;
nSubbands = 4;
carrierFrequency = centerFrequency - bandwidth * (1 - 1 / nSubbands) / 2: bandwidth / nSubbands: centerFrequency + bandwidth * (1 - 1 / nSubbands) / 2;
distance = 10;
pathlossExponent = 2;
nRealizations = 1e3;
fadingType = 'selective';
%% Harvester
% assumptions: antenna impedance = 50 ohms, ideality factor = 1, thermal voltage = 25.85 mV
beta2 = 9.6712e2;
beta4 = 6.0304e6;
tolerance = 5e-3;
%% Variables
Variable.nUsers = [2 3];
Variable.nCandidates = [1 5 50 500];

%% Transceiver
eirp = db2pow(36 - 30);
nCandidates = 50;
nTxs = 20;
txPower = eirp / nTxs;
%% Channel
centerFrequency = 2.4e9;
bandwidth = 1e7;
nSubbands = 8;
[carrierFrequency] = carrier_frequency(centerFrequency, bandwidth);
distance = 10;

nRealizations = 2;
fadingType = 'selective';
%% Harvester
% assumptions: antenna impedance = 50 ohms, ideality factor = 1, thermal voltage = 25.85 mV
beta2 = 9.6712e2;
beta4 = 6.0304e6;
tolerance = 1e-3;
%% Variables
Variable.nUsers = [2 5];
%% Transceiver
eirp = db2pow(36 - 30);
nCandidates = 50;
%% Channel
centerFrequency = 2.4e9;
bandwidth = 1e7;
distance = 10;
pathlossExponent = 2;
nRealizations = 1;
fadingType = 'selective';
%% Harvester
% assumptions: antenna impedance = 50 ohms, ideality factor = 1, thermal voltage = 25.85 mV
beta2 = 9.6712e2;
beta4 = 6.0304e6;
tolerance = 1e-3;
%% Variables
Variable.nTxs = [4, 20];
Variable.nSubbands = [4, 8];
Variable.nUsers = 2 : 5;

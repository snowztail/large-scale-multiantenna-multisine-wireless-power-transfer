%% User
nUsers = 1;
%% Transceiver
powerBudget = 0.5;
%% Channel
centerFrequency = 2.4e9;
bandwidth = 1e7;
pathlossExponent = 2;
nRealizations = 2e2;
fadingType = 'selective';
%% Harvester
% assumptions: antenna impedance = 50 ohms, ideality factor = 1, thermal voltage = 25.85 mV
beta2 = 9.6712e2;
beta4 = 6.0304e6;
tolerance = 1e-3;
%% Variables
Variable.distance = 10 : 2 : 20;
Variable.nTxs = [8 16 16 32 32];
Variable.nSubbands = [8 8 16 16 32];

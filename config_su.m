clc; clear; close all;
%% Transceiver

%% Channel
% center frequency
centerFrequency = 2.4e9;
% bandwidth
bandwidth = 1e7;
% number of frequency bands
subband = 8;
%
fs = 20e6; % Channel model sampling frequency equals the channel bandwidth
% tgnChan = wlanTGnChannel('CarrierFrequency', centerFrequency, 'SampleRate', bandwidth, 'LargeScaleFadingEffect', 'Pathloss', 'DelayProfile', 'Model-E');
% Create and configure the channel
tgnChannel = wlanTGnChannel;
tgnChannel.DelayProfile = 'Model-B';
tgnChannel.NumTransmitAntennas = 2;
tgnChannel.NumReceiveAntennas = 2;
tgnChannel.TransmitReceiveDistance = 10; % Distance in meters for NLOS
tgnChannel.LargeScaleFadingEffect = 'None';

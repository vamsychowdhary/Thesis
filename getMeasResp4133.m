function [mic_4133_calib_dB,mic_4133_ph,fn] = getMeasResp4133

% This function output the measured and calibrated magnitude and phase
% responses of B&K 4133 Microphone. The data is take from 31220 course at
% DTU.

load('part1.mat','H_21_ref','fn');

load('part3final.mat','specn3b2_4133','specn3c_4133');

% Measure the Noise
mic_4133_noise = (specn3b2_4133(:,2)./specn3b2_4133(:,1))./H_21_ref;

% Remove the Noise from Signal
mic_4133 = (specn3c_4133(:,2)./specn3c_4133(:,1))./H_21_ref - mic_4133_noise;

% Calculate the calibrated frequency response
sensitivity = 39.95e-3;
M = 12.6e-3;
mic_4133_calib = (mic_4133./sensitivity).*M;
mic_4133_calib_dB = 20*log10(abs(mic_4133_calib));
mic_4133_ph = atan2d(imag(mic_4133_calib),real(mic_4133_calib));

end
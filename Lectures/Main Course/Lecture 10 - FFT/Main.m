%% LECTURE 10 - Convolution / FFT Methods - Noé Debrois - 30/10/2024
% This code implements the calling of the CONV function to price a call 
% option of type up and out with a certain barrier.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Model parameters :
param.rf = 0.05; % Risk-free rate
param.q = 0.02;     % Dividend
param.distr = 1; % Normal distribution (if 2 : NIG)
param.m = 0.0;     % Drift
param.s = 0.2;   % Volatility

% Option parameters :
S_0 = 100;         % Spot price
K = 95;           % Strike
Barrier = 80;   % Barrier value
N = 2^14;        % Grid of log-returns (log(S_T / S_0))

% Time parameters :
param.T = 1;                % Maturity
Ndate = 12;                 % Number of monitorings till maturity
param.dt = param.T / Ndate; % Monthly monitoring


[S,v] = CONV(S_0, K, Ndate, N, Barrier, param);
price = interp1(S, v, S_0, 'spline') % S_0, our spot price, is not 
% necessarily contained exactly in the S vector returned by CONV fct. So we
% need to interpolate it.
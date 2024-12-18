clear; close all;
addpath("CHAR_EXP/");

% Test parameters
S_0 = 100;       % Initial spot price
K = 95;          % Strike price
Ndate = 12;      % Number of time steps (dates)
N = 512;         % Grid size for Fourier transform
Barrier = 80;    % Barrier level (for Down&Out option)
param.sigma = 0.2;  % Volatility
param.dt = 1/252;   % Time step (daily)
param.rf = 0.05;    % Risk-free rate
param.q = 0.02;     % Dividend yield
param.distr = 1;    % Black-Scholes (Normal distribution)

% MERTON (2)
param.mu = 0.01;
param.delta = 0.2;
param.lambdaK = 2;

% KOU (3)
param.p = 0.6;
param.lambdap = 10; 
param.lambdam = 3; 

% NIG (4)
param.kNIG = 1; % kNIG = 0.001 -> Like B&S.
param.theta = 0.05;

% VG (5)
param.kVG = 0.05; % Variance of the subordinator

[S,v] = CONV_CHAR_FCT(S_0, K, Ndate, N, Barrier, param);
price = interp1(S, v, S_0, 'spline')
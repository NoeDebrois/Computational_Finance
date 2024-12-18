addpath("CHAR_EXP/");

% Test parameters
S_0 = 100;       % Initial spot price
K = 95;          % Strike price
Ndate = 100;     % Number of time steps (dates)
N = 512;         % Grid size for Fourier transform
Barrier = 80;    % Barrier level (for Down&Out option)
param.sigma = 0.2;  % Volatility
param.dt = 1/252;   % Time step (daily)
param.rf = 0.05;    % Risk-free rate
param.q = 0.02;     % Dividend yield
param.distr = 1;    % Black-Scholes (Normal distribution)

% Run the function
[S, v] = CONV_CHAR_FCT(S_0, K, Ndate, N, Barrier, param);

% Display results
disp('Spot Prices and Corresponding Option Prices:');
disp([S, v]);

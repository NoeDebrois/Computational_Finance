%% Pricing script of a CALL BARRIER OPTION with DISCRETE MONITORING:
clear; close all;
addpath("CHAR_EXP/");
%-------------------------------------------------------------------------%
%% Call option parameters:
S_0 = 100;       % Initial spot price
K = 95;          % Strike price
param.T = 1;     % Maturity

param.rf = 0.05; % Risk-free rate
param.q = 0.02;  % Dividend yield

Ndate = 120;      % Number of monitoring dates (12: Monthly)
Barrier = 80;    % Barrier level (for Down&Out option)
param.dt = param.T / Ndate; % Monitoring (1/12: monthly)

%% Simulation parameters:
param.distr = 5;    % Model choice :
                    % B&S: 1 ; Merton: 2 ; Kou: 3 ; NIG: 4 ; VG: 5 ;

N = 2^14;           % Grid size for Fourier transform

%% Model parameters: choose your parameters depending on your model
if param.distr == 1 % B&S
    % sigma = volatility of the model
    param.sigma = 0.2;

elseif param.distr == 2 % Merton
    % sigma   = BM volatility (as in GBM)
    % mu      = jump drift (meansize of jumps)
    % delta   = jump vol (std of jumps)
    % lambdaK = jump time intensity
    param.sigma   = 0.2;
    param.mu      = 0.01;
    param.delta   = 0.2;
    param.lambdaK = 2;

elseif param.distr == 3 % Kou
    % sigma   = BM volatility (as in GBM)
    % p       = prob of positive jump
    % lambdap = pos jump intensity
    % lambdam = neg jump intensity
    % lambdaK = jump time intensity
    param.sigma   = 0.2;
    param.p       = 0.6;
    param.lambdap = 10;
    param.lambdam = 3;
    param.lambdaK = 2;

elseif param.distr == 4 % NIG
    % sigma = volatility of the BM
    % theta = drift of the BM
    % kNIG  = variance of the subordinator
    param.sigma = 0.2;
    param.theta = 0.05;
    param.kNIG  = 1; % kNIG->0^+ : like B&S.

elseif param.distr == 5 % VG
    % sigma = volatility of the BM
    % theta = drift of the BM
    % k     = variance of the subordinator
    param.sigma = 0.2;
    param.theta = 0.05;
    param.kVG   = 0.05;

%elseif param.distr == 6
%elseif param.distr == 7
end

[S,v] = CONV_CHAR_FCT(S_0, K, Ndate, N, Barrier, param);
price = interp1(S, v, S_0, 'spline');
disp(price);
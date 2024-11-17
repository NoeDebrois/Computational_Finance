%% LECTURE 16 - Heston Model - Noé Debrois - 16/11/2024
% This code implements the simulation the Heston model via the Euler scheme
% (not the article from Andersen). So it is like in CIR_Simulation.m, but
% with also the simulation of X (using V). Of course, V can be negative,
% even when Feller condition is satisfied, which is the main problem of
% this method. We need more advanced schemes to solve this (cf Andersen).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
%% Inputs :
% Variance : 
V0 = 0.1; % Variance at time 0.
theta = 0.3; % Mean.
k = 0.1; % Speed of mean reversion.
epsilon = 0.2; % Vol-of-var.

% Price :
r = 0.02;
S0 = 100; 

% Correlation between price and variance :
rho = -0.6;

% Calibration :
% -> r, S0 : observed on the market ;
% -> [V0, theta, k, epsilon, rho] have to be calibrated (from European Call
% exploiting Carr-Madan --> characteristic function of Heston known).

Feller_condition = (epsilon^2 - 2 * k * theta)
if Feller_condition <= 0
    disp('Feller Condition is satisfied')
else
    disp('Feller Condition is NOT satisfied')
end

% Discretization :
Nsim = 1e4; T = 1; Nsteps = 25;

% MC simulation :
dt = T / Nsteps; t = linspace(0, T, Nsteps + 1); % Time grid.
x = zeros(Nsim, Nsteps + 1); x(:, 1) = log(S0); % Log price process.
V = zeros(Nsim, Nsteps + 1); V(:, 1) = V0; % Variance process.
%
%% Mean & variance-covariance matrix of normal distribution to be sampled :
mu = [0; 0]; VC = [1 rho; rho 1]; % Covariance matrix between x and V.
for i=1:Nsteps
    Z = mvnrnd(mu, VC, Nsim); % Multivariate normal random numbers.
    x(:, i + 1) = x(:, i) + (r - V(:, i) / 2) * dt +...
        sqrt(V(:, i) * dt) .* Z(:, 1);
    %
    V(:, i + 1) = V(:, i) + k * (theta - V(:, i)) * dt +...
        epsilon * sqrt(V(:, i) * dt) .* Z(:, 2);
end

% Computation of the price from the log-price :
S = exp(x);
%
%% Plot :
figure; hold on
subplot(1, 2, 1); plot(t, S); title('Price');
subplot(1, 2, 2); plot(t, V); title('Variance');
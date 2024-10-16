%% Lecture 8 - Monte Carlo Simulation - Noé Debrois - 16/10/2024
% MC, Variance Gamma (VG), Variance reduction technique : Control Variable (CV).
% Pricing an European Call under VG by MC with CV.
clear
clc
%
%% Parameters :
% MC parameter :
Nsim = 1e6; Nsim2 = 1000;
% Market parameters :
S0 = 100; r = 0.02;
% Contract parameters: 
K = 105; T = 2;
% Model (VG) parameters :
sigma = 0.6; theta = 0.2; k = 0.5;
%
%% 1. Simulation :
% For details, cf course on Control Variable MC (page 83).
%

% Here we choose f = ST = S0*exp((r-sigma^2/2)T+sigmasqrt(T)Z) (cf notes).
% So, under the RISK-NEUTRAL measure : E[f(Z)] = S0*exp(rT) :
Ef = S0 * exp(r * T);

% Characteristic exponent for VG (cf VG_NIG.pdf) :
char_exp = @(u) -log(1 + u.^2 * sigma^2 * k / 2 - 1i * theta * k * u) / k;
% Drift :
drift = r - char_exp(-1i);

% A. Small simulation (Nsim2) to compute alpha :
%
% cf notes (page 83).
%
% - We don't know exactly alpha, so we have to simulate it. In order not to
% be too computationally-intense, we'll do a small number of simulations.
dS = k * icdf('gamma', rand(Nsim2, 1), T / k, 1);
XT = 0 + drift * T + theta * dS + sigma * sqrt(dS) .* randn(Nsim2, 1);
f = S0 * exp(XT); % f(Z) = S0*exp(XT) = S0*exp((r-sigma^2/2)T+sigmasqrt(T)Z)
g = exp(-r * T) * max(f - K, 0); % g(Z) = exp(-rT)*max(f(Z)-K, 0)
A = cov(g, f);
alpha = - A(1,2) / A(2,2); % alpha = - cov(f(Z), g(Z)) / Var(f(Z))

% B. Compute the price (now that we have a (bad) estimation of alpha) :
dS = k * icdf('gamma', rand(Nsim, 1), T / k, 1);
XT = 0 + drift * T + theta * dS + sigma * sqrt(dS) .* randn(Nsim, 1);
f = S0 * exp(XT); 
g = exp(-r * T) * max(f - K, 0);
%
%% 2. Price and Confidence Interval :
% In CV we study g(Z) + alpha * (f(Z) - g(Z)) :  
[price, ~, CI] = normfit(g + alpha * (f - Ef))
% Let's try Carr Madan method to compare :
priceCM = FFT_CM_Call_VG(K, [sigma,theta,k], T, r, S0)
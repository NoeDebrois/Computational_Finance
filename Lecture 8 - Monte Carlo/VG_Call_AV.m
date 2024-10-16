%% Lecture 8 - Monte Carlo Simulation - Noé Debrois - 16/10/2024
% MC with variance reduction : Antithetic Variable (AV), Variance Gamma (VG).
% Pricing an European Call under VG by MC with AV.
clear
clc
%
%% Parameters :
% MC parameter :
Nsim = 1e6;
% Market parameters :
S0 = 100; r = 0.02;
% Contract parameters: 
K = 105; T = 2;
% Model (VG) parameters :
sigma = 0.6; theta = 0.2; k = 0.5;
%
%% 1. Simulation :
% For details on VG, see VG_NIG.pdf ;
% For details on how to simulate VG, see Simulate_VG_NIG.pdf :
% "ALGORITHM 6.11 : Simulating a VG process on a fixed time grid" (p.4).
%

% Characteristic exponent for VG (cf VG_NIG.pdf) :
char_exp = @(u) -log(1 + u.^2 * sigma^2 * k / 2 - 1i * theta * k * u) / k;
% Drift :
drift = r - char_exp(-1i);

% STEP 1 : simulate 1 independent gamma variable with parameter dt.
% Here M = 1 since it's not path dependant : we only need the final price !
%
% dS = DeltaS ~ gamma variable with parameter dt = T.
% Here dt = T because dt = T / M = T / 1 since it's not path dependant !
%
dS = k * icdf('gamma', rand(Nsim, 1), T / k, 1);

% STEP 2 : simulate M (= 1) RV ~ N(0, 1).
Z = randn(Nsim, 1); % size : Nsim x 1, for the MC.
%
% And compute :
% X(t{i+1}) = X(t{i}) + drift_part 
%             + theta * DeltaS 
%             + sigma * sqrt(DeltaS) * Z
%
XT = 0 + drift * T + theta * dS + sigma * sqrt(dS) .* Z;
% Once again, here M = 1 : this option is not path dependant so we only
% need the final price.

% Build ST1 and ST2 using the Antithetic Variable (AV) technique :
ST1 = S0 * exp(XT);
% AV :
XT = 0 + drift * T + theta * dS - sigma * sqrt(dS) .* Z;
ST2 = S0 * exp(XT);

% Check if the prices are simulated correctly :
%
% Theoretical expectation: under the risk-neutral measure, the expected 
% value of the asset price at time T is given by : S0*exp(rT), where 
% - r : the risk-free rate ;
% - T : the time horizon.
% This comes from the assumption that the stock price follows the 
% discounted martingale property under the risk-neutral measure.
%
% So, in average, ST - SO*exp(rT) should be equal to zero.
%
% Let's try for both ST1 and ST2 :
[check1, ~, CI_check] = normfit(ST1 - S0 * exp(r * T))
[check2, ~, CI_check] = normfit(ST2 - S0 * exp(r * T))
%
%% 2. Discounted Payoff :
% EU Call Option, with both simulated prices (AV technique) :
discpayoff1 = exp(-r * T) * max(ST1 - K, 0);
discpayoff2 = exp(-r * T) * max(ST2 - K, 0);
%
%% 3. Price and Confidence Interval :
% EU Call Option, we average the two discounted payoffs (AV technique) :
[price, ~, CI] = normfit(1/2 *(discpayoff1 + discpayoff2))
% And Carr Madan method to compare :
priceCM = FFT_CM_Call_VG(K, [sigma,theta,k], T, r, S0)
% Normally, you should see a better convergence here than w/ VG_Call_AV.m.
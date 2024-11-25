%% Lecture 8 - Monte Carlo Simulation - Noé Debrois - 16/10/2024
% MC, Variance Gamma (VG), 2.
% Pricing Knock&Out Call Option under VG by MC.
clear
clc
%
%% Parameters :
% MC parameter :
Nsim = 1e6;
% Market parameters :
S0 = 100; r = 0.02;
% Contract parameters :
K = 105; T = 2; U = 130; L = 80;
% Fixed time grid monitoring :
M = round(12 * T); % [monthly -> 12 ; weekly -> 52]
% Model (VG) parameters :
sigma = 0.6; theta = 0.2; k = 0.5;
%
%% 1. Simulation :
% For details on VG, see VG_NIG.pdf ;
% For details on how to simulate VG, see Simulate_VG_NIG.pdf :
% "ALGORITHM 6.11 : Simulating a VG process on a fixed time grid" (p.4).
%
dt = T / M; % Time step
X = zeros(Nsim, M + 1); % Initialisation of the process X

% Characteristic exponent for VG (cf VG_NIG.pdf) :
char_exp = @(u) -log(1 + u.^2 * sigma^2 * k / 2 - 1i * theta * k * u) / k;
% Drift :
drift = r - char_exp(-1i);

% STEP 1 : simulate M independent gamma variables with parameter dt.
%
% dS = DeltaS [1 to M] ~ independant gamma variables with parameter dt
% (because dt = [t_{i}-t_{i-1}] / k) ; and multiply by k :
%
dS = k * icdf('gamma', rand(Nsim, M), dt / k, 1); % size : Nsim x M. 

% STEP 2 : simulate M i.i.d ~ N(0, 1) RVs.
Z = randn(Nsim, M); % size : Nsim x M, for the MC.
%
% And compute :
% X(t{i+1}) = X(t{i}) + drift_part 
%             + theta * DeltaS(t{i}) 
%             + sigma * sqrt(DeltaS(t{i})) * Z(t{i})
%
for j=1:M
    X(:,j+1) = X(:,j) + drift * dt + theta * dS(:,j) + sigma * sqrt(dS(:,j)) .* Z(:,j);
end

% Build S :
S = S0 * exp(X);
clear X dS Z;

% Check if the price is simulated correctly :
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
disp('Risk-Neutrality check: in average, ST - SO*exp(rT) must be zero !')
[check, ~, CI_check] = normfit(S(:,end) - S0 * exp(r * T))
%
%% 2. Discounted Payoff :
% Knock&Out Call Option :
discpayoff = exp(-r * T) * max(S(:,end) - K, 0) .* (min(S,[],2) > L) .* (max(S,[],2) < U);
% (min(S,[],2) > L) = 1 if it is true that, for the given MC path, the
% price was always above the lower bound. Analogous for the upper bound.
%
%% 3. Price and Confidence Interval :
% Knock&Out Call Option :
[price, ~, CI] = normfit(discpayoff) % MC CI
%% LECTURE 16 - Heston Model - Noé Debrois - 22/11/2024
% This code implements the QE Scheme for Heston model, from the article of
% Andersen : "Efficient Simulation of the Heston Stochastic Volatility 
% Model". -> A fully commented function version is available : Heston_QE.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%--- PARAMETERS -----------------------------------------------------------
V0 = 0.01; % Variance at time 0.
k = 2; % Speed of mean reversion.
theta = 0.002; % Long-term mean.
epsilon = 0.2; % Vol-of-var.

Nsim = 10; % Nb of sims.
M = 100; T = 1; dt = T / M; % Time grid.

V = zeros(Nsim, M + 1); V(:, 1) = V0; % Variance.

% Choose Psi_cutoff = 1.5 as suggested in article :
Psi_cutoff = 1.5;

%--- DISCRETIZE V --------------------------------------------------------
for i=1:M
    % STEP 1 & 2
    % Calculate m,s,Psi
    m = theta + (V(:, i) - theta) * exp(-k * dt);
    m2 = m.^2;

    s2 = V(:, i) * epsilon^2 * exp(-k * dt) * ...
        (1 - exp(-k * dt)) / k + ...
        theta * epsilon^2 * ...
        (1 - exp(-k * dt))^2 / (2 * k);
    s = sqrt(s2);
    
    Psi = (s2)./(m2);
    
    % STEP 3,4,5
    % Depending on Psi, use exp or quad scheme to calculate next V
    index = find(Psi > Psi_cutoff);% Detection for E scheme.

    % Exponential approximation, used for Psi>Psi_cutoff)
    % PDF of V(t+dt) is p*delta(0) + (1-p) * (1-exp(-beta x )
    % thus a prob. mass in 0 & an exp tail after that
    p_exp = (Psi(index) - 1) ./ (Psi(index) + 1 ); % Prob. mass in 0, Eq 29
    beta_exp = (1 - p_exp) ./ m(index);	% Exponent of exp. density tail, Eq 30
    % Gets x from inverse CDF applied to uniform U
    U = rand(size(index));
    V(index, i + 1) = (log((1 - p_exp) ./ (1 - U)) ./ beta_exp) .* (U > p_exp);

    index = find(Psi <= Psi_cutoff); % Detection for Q scheme.
    % Quadratic approx, used for 0<Psi<Psi_cutoff
    % V(t+dt) = a(b+Zv)^2, Zv~N(0,1)
    invPsi = 1 ./ Psi(index);
    % Eq 27 :
    b2_quad = 2 * invPsi - 1 + sqrt(2 * invPsi) .* sqrt(2 * invPsi - 1);
    % Eq 28 :
    a_quad = m(index) ./ (1 + b2_quad);	
    % Eq.23 :
    V(index,i+1) = a_quad .* (sqrt(b2_quad) + randn(size(index))) .^2; 
end

%--- PLOT -----------------------------------------------------------------
t = linspace(0, T, M + 1);
plot(t, V)
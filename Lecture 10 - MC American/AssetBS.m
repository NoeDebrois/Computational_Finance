%% LECTURE 10 - MC for American Options - Noé Debrois - 28/10/2024
% Function to compute the asset price for the MC simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function S = AssetBS(r, sigma, S0, T, M, Nsim)
% Simulates an asset price according to the B&S Framework.
X = zeros(Nsim,M+1);
dt = T / M;
for i=1:M
    X(:, i + 1) = X(:, i) + (r - sigma^2 / 2) * dt +...
        sigma * sqrt(dt) * randn(Nsim, 1);
end
S = S0 * exp(X);
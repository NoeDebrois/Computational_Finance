%% Simulation of the Merton model V2
% This time we simulate the jump at each time step (in the V1 we simulated
% all the jumptimes before the loop).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
% General parameters :
S0 = 1; % Initial value
mu = 0.05; % Drift
sigma = 0.4; % Volatility
lambda = 2; % Poisson intensity/rate
% Jumpsize parameters : Jumpsize ~ Normal(muJ, deltaJ^2) :
muJ = 0.01; % mean of the size of the jumps
deltaJ = 0.2; % std deviation of the size of the jumps
% Simulation parameters :
T = 2; % Maturity
M = 100; % Number of steps in time
dt = T/M; % Time step
% 
%% Simulation :
% Initialisations of X (in the exp) & Z (vector of N(0, 1) RV) :
X = zeros(M+1,1);
Z = randn(M,1);
%
% cf "Simulate_Jump_Diffusion.pdf" : ALGORITHM 6.1.
for i=1:M
    % 1) Simulate the continuous part of the path :
    %
    X(i+1) = X(i) + mu * dt + sigma * sqrt(dt) * Z(i);
    %
    % 2) Add jumps :
    %
    Ndt = poissrnd(lambda * dt);
    % Ndt represents the number of jumps in the current time step.
    % It is simulated from a Poisson distribution with parameter lambda * dt, 
    % which represents the intensity of jumps occurring during this time interval.
    if Ndt == 0
        J = 0; % No jump occurs during this time step, so J is set to 0.
    else
        J = sum(muJ + deltaJ * randn(Ndt,1));
        % If Ndt > 0, we have one or more jumps during this time step.
        % Each jump is modeled as a normally distributed random variable 
        % with mean muJ and standard deviation deltaJ.
        % We sum up the effects of all the jumps that occurred in this time step.
    end
    X(i+1) = X(i+1) + J; % We add the simulated jump size for this time step.
end
%% Plot :
figure;
plot(S0 * exp(X));
title("One simulated path via Merton model, V2 :")
    
    
    
    
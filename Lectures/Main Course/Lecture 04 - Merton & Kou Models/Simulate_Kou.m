%% Simulation of the Kou model :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
% S_T = S0 * exp(X_T) with :
% (X_t)_t ~ Kou model : Jumpsize ~ Asymmetric Exp(param) dist (param =
% lambda+ for positive jumps, param = lambda- for negative jumps, p is the
% probability of positive jump).
%
% General parameters :
S0 = 1; % Initial value
mu = 0.05; % Drift
sigma = 0.4; % Volatility
lambda = 2; % Poisson intensity
% Jumpsize parameters : Jumpsize ~ Exp(BLABLA) :
p = 0.6; % Probability of having a positive jump
lambda_plus = 10; % Parameter of Exp() for POSITIVE jumps
lambda_minus = 3; % Parameter of Exp() for NEGATIVE jumps
% WARNING : the mean size of jump is 1/lambda_whatever !! so don't take
% small values of lambda_whatever...
%
% Simulation parameters :
T = 2; % Maturity
M = 100; % Number of steps in time
dt = T/M;
%% Simulation :
% Initialisations of X (in the exp) & Z (vector of N(0, 1) RV) :
X = zeros(M+1,1);
Z = randn(M,1);
% 
% CONDITIONAL SIMULATION of jump times (the easiest way) :
% cf "Simulate_Jump_Diffusion.pdf" : ALGORITHM 6.2
% - Generates random nb from Poisson distribution of parameter lambda*t:
NT = poissrnd(lambda*T);
% - Generates and orders jump times chronologically :
jumpT = sort(rand(1,NT)*T); % sort(a vector of U() RV of length NT, rescaled between 0 and T)
% 
for i=1:M
    % 1) Simulate the continuous part of the path :
    %
    X(i+1) = X(i) + mu * dt + sigma * sqrt(dt) * Z(i);
    % 2) Add jumps : we use CONDITIONAL SIMULATION of jump times (cf above)
    % - Loop on all the jump times ;
    % - Check if there are jumps in ((i-1) * dt, i * dt] ;
    % - If yes, add the simulated jump. -> Kou simulated jump.
    for j=1:NT % Loop on all the jump times
        % Does the jump occur between the 2 timesteps we consider ?
        if jumpT(j) > (i-1) * dt && jumpT(j) <= i * dt
            %
            % Yes, there is a jump : let's simulate its sign + size :
            %
            u = rand; % RV ~ U((0,1)) to know which type of jump (+ or -).
            if u < p 
               %
               % Positive jump :
               %
               jumpSize = exprnd(1/lambda_plus);% jumpsize ~ Exp(lambda_+)
            else
               %
               % Negative jump :
               %
               jumpSize = -exprnd(1/lambda_minus);% jumpsize ~ Exp(lambda_-)
            end  
            % We add the simulated jump size :
            X(i+1) = X(i+1) + jumpSize;
        end
    end
end
%% Plot :
figure;
plot(S0 * exp(X));
title("One simulated path via Kou model :")
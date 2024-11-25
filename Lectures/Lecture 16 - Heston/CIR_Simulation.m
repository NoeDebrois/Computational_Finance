%% LECTURE 16 - Heston Model - Noé Debrois - 16/11/2024
% This code implements the sampling of CIR using Euler Scheme. We can see,
% if we zoom, non-zero probability to go below 0, whereas the theory says
% it's impossible : NUMERICAL ERROR (see the warning in the command line
% about "Imaginary Part". So we prevent it by a max(y, 0)...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
% t_i = i * dt, let's discretize the SDE from (y_t)_t :
% y_{i+1} = y_{i} + lambda * (eta - y_{i}) * dt + beta * sqrt(y_{i}) *
% sqrt(t) * Z, Z ~ N(0, 1) [EULER SCHEME].
%
%% Parameters :
y0 = 0.02;
beta = 0.3;
lambda = 1;
eta = 0.05;

if beta^2 - 2 * lambda * eta <= 0
    disp('Feller condition satisfied')
    % P(y_t = 0) = 0 & P(y_t < 0) = 0
else
    disp('Feller condition NOT satisfied')
    % P(y_t = 0) > 0 & P(y_t < 0) = 0
end

T = 1;
M = 100;
dt = T / M;
Nsim = 50; 
%
%% Sampling CIR using Euler Scheme :
y = zeros(Nsim, M + 1); 
y(:, 1) = y0; % Initial condition

for i=1:M
    y(:, i + 1) = y(:, i) + lambda * (eta - y(:, i)) * dt +...
        beta * sqrt(y(:, i)) * sqrt(dt) .* randn(Nsim, 1);
    y(:, i + 1) = max(y(:, i + 1), 0); % Prevent from falling below 0 !
end
%
%% Plot :
t = linspace(0, T, M + 1);
plot(t, y);
title("Realisations of (y_t)_t (CIR) using Euler Scheme")
%% LECTURE 12 - Theta Method - Noé Debrois - 04/11/2024
% This code implements the pricing of a Plain Vanilla Call Option, under
% B&S model, using the Theta Method on the logprice-PDE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
%% 1. Parameters & Grids :
% Model & Option parameters :
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Method parameter :
theta = 0.5; 
% 0 : Implicit Euler, 0.5 : Crank-Nicholson, 1 : Explicit Euler

% Space grid :
N = 2000; % Number of points in the space dimension.
Smin = 0.1 * S0; Smax = 3 * S0; % Found by experience : we assume that the
% spot price has a negligible probability of falling below 10% of S0, and
% of rising above 3 * S0.
xmin = log(Smin / S0); xmax = log(Smax / S0); % x lives in [xmin, xmax]=
% [log(Smin/S0),log(Smax/S0)].
dx = (xmax - xmin) / N; % Space grid step.
x = xmin + (0:N) * dx; % From x0 = xmin to xN = xmax (in total : N+1 pts).
% x = linspace(xmin, xmax, N+1);

% Time grid :
M = 1000;       % Number of points in the time dimension.
dt = T / M;     % Time grid step.
t = (0:M) * dt; % t lives in [0, T], divided in M subintervals : M+1 pts.
%
%% 2. Matrix Construction :
% We write the theta method in matrix form, cf my notebook about PDE.

% [A, B, C] :
A = (1 - theta) * dt * (-(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
B = -1 + dt * (1 - theta) * (-sigma^2 / (dx^2) - r);
C = dt * (1 - theta) * ((r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
%
Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2); % M1 in my notes.
Mat(1, 1) = 1;
for i=2:N
    Mat(i, [i - 1, i, i + 1]) = [A B C];
end
Mat(end, end) = 1;

% [Atilde, Btilde, Ctilde] :
Ah = -(theta) * dt * (-(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
Bh = -1 -dt * (theta) * (-sigma^2 / (dx^2) - r);
Ch = -dt * (theta) * ((r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
%
Mat_rhs = spalloc(N + 1, N + 1, 3 * (N - 1)); % M2 in my notes.
for i=2:N
    Mat_rhs(i, [i - 1, i, i + 1]) = [Ah Bh Ch];
end
%
%% 3. Backward in time loop :
c = max(S0 * exp(x') - K, 0); % @ Maturity : "v(T, x) = (S0 * e^x - K)^+".

rhs = zeros(N + 1, 1); % BCj in my notes.

for j=M:-1:1 % BACKWARD from t_M to t_0.
    % We know c t_j -> we compute cnew t_{j-1}.
    rhs = Mat_rhs * c;
    rhs(1) = 0;
    rhs(end) = Smax - K * exp(-r * (T - (j - 1) * dt)); % BC at tj, for x=xmax
    c = Mat \ rhs;
end
%
%% Plot "call price at t=0" VS "S at t=0" (for all S in [Smin, Smax]) :
figure
plot(S0 * exp(x), c);
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0');
%
%% Compute the price at t = 0 :
price = interp1(x, c, 0, 'spline') % Interpolation of c at x = 0, i.e S = S0
% (S0 being the one we the one we initialized in "Model parameters").
price_ex = blsprice(S0, K, r, T, sigma) % For comparison : B&S price.
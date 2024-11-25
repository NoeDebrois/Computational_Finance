%% LECTURE 11 - Theta Method for EU Call - Noé Debrois - 29/10/2024
% Script to price a EU Call Option, under B&S model, using PRICE PDE.
% WARNING : NOT log-price PDE ! cf Lecture notes, page 23 PDE.pdf.
% We implement Finite Difference Method using THETA METHOD SCHEME.
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
dS = (Smax - Smin) / N; % Space grid (price grid) step.
S = Smin + (0:N) * dS; % Price grid.

% Time grid :
M = 1000;       % Number of points in the time dimension.
dt = T / M;     % Time grid step.
t = (0:M) * dt; % t lives in [0, T], divided in M subintervals : M+1 pts.
%
%% 2. Matrix Construction :
% We write the theta method in matrix form, cf my notebook about PDE.


% WARNING : with the PRICE PDE (contrary to logprice PDE), A, B & C depend
% on i. Hence the following functions to build A, B, C for each S_i :
A = @(S) (1 - theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
B = @(S) -1 + dt * (1 - theta) * (-sigma^2 * S^2 / (dS^2) - r);
C = @(S) dt * (1 - theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
% cf page 26 of notes PDE.pdf.

Mat = spalloc(N + 1, N + 1, 3 * (N - 1) + 2);
Mat(1, 1) = 1;
for i=2:N
    Mat(i, [i - 1, i, i + 1]) = [A(S(i)) B(S(i)) C(S(i))];
end
Mat(end, end) = 1;

% [Atilde, Btilde, Ctilde] (once again, with PRICE PDE, they depend on i) :
Ah = @(S) -(theta) * dt * (-r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));
Bh = @(S) -1 - dt * (theta) * (-sigma^2 * S^2 / (dS^2) - r);
Ch = @(S) -dt * (theta) * (r * S / (2 * dS) + sigma^2 * S^2 / (2 * dS^2));

Mat_rhs = spalloc(N + 1, N + 1, 3 * (N - 1));
for i=2:N
    Mat_rhs(i, [i - 1, i, i + 1])=[Ah(S(i)) Bh(S(i)) Ch(S(i))];
end
%
%% 3. Backward in time loop :
c = max(S' - K, 0); % @ Maturity : "v(T, x) = (ST - K)^+". 

rhs = zeros(N + 1, 1); % BCj in my notes.

for j=M:-1:1 % BACKWARD from t_M to t_0.
    % We know c t_j -> we compute cnew t_{j-1}.
    rhs = Mat_rhs * c;
    rhs(1) = 0;
    rhs(end) = Smax - K * exp(-r * (T - (j - 1) * dt)); % BC at tj, for S=Smax
    c = Mat \ rhs;
end
%
%% Plot "call price at t=0" VS "S at t=0" (for all S in [Smin, Smax]) :
figure
plot(S, c); 
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0')
%
%% Compute the price at t = 0 :
price = interp1(S, c, S0, 'spline') % Interpolation of c at spot price S0 (fixed above).
price_ex = blsprice(S0, K, r, T, sigma) % For comparison : B&S price.
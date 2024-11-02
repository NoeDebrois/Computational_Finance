%% LECTURE 11 - Implicit Euler for EU Call - Noé Debrois - 29/10/2024
% Script to price a EU Call Option, under B&S model, using log-price PDE.
% We implement Finite Difference Method using IMPLICIT EULER SCHEME.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
%
%% 1. Parameters & Grids :
% Model parameters :
S0 = 1; K = 0.95; r = 0.001; sigma = 0.5; T = 1;

% Space grid :
N = 2000; % Number of points in the space dimension.
Smin = 0.1 * S0; % Found by experience : we assume that the spot price has
Smax = 3 * S0;   % a negligible probability of falling below 10% of S0, and
                 % of rising above 3 * S0.
xmin = log(Smin/S0); % x lives in [xmin, xmax]=[log(Smin/S0),log(Smax/S0)].
xmax = log(Smax/S0);
dx = (xmax - xmin) / N; % Space grid step.
x = xmin + (0:N) * dx;  % From x0 = xmin to xN = xmax (in total : N+1 pts).
% x = linspace(xmin, xmax, N+1);

% Time grid :
M = 1000;   % Number of points in the time dimension.
dt = T / M; % Time grid step.
t = (0:M) * dt; % t lives in [0, T], divided in M subintervals : M+1 pts.
%
%% 2. Implicit Euler - Matrix Construction :
% Implicit Euler Scheme : by moving V_{j+1, i} (KNOWN) on the right :
% V_{j, 0} = 0 [BOUNDARY CONDITION] ;
% A * V_{j, 0} + B * V_{j, 1} + C * V_{j, 2} = (-1/dt) * V_{j+1, 1} ;
% [...]
% A * V_{j, N-2} + B * V_{j, N-1} + C * V_{j, N} = (-1/dt) * V_{j+1, N-1} ;
% V_{j, N} = Smax - K * exp(-r * (T - j * dt)) [BOUNDARY CONDITION].
%
% Where V_{j, i-1}, V_{j, i} & V_{j, i+1} are UNKNOWN.
%
% And :
%       A = - (r - sigma^2 / 2) / (2 * dx) + (sigma^2 / 2) / dx^2 ;
%       B = - 1 / dt - sigma^2 / dx^2 - r ;
%       C = (r - sigma^2 / 2) / (2 * dx) + (sigma^2 / 2) / dx^2.
%
A = dt * (-(r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));
B = -1 + dt * (- sigma^2 / (dx^2) - r);
C = dt * ((r - sigma^2 / 2) / (2 * dx) + sigma^2 / (2 * dx^2));

% Because each equation has 3 unknown : it's not enough to come up with one
% solution. But with the boundary conditions we can do it. This corresponds
% to : "Mat * V = Rhs", with :
% - First row : [1 0 ... 0] ;
% - Second row : [A B C 0 ... 0] ;
% - Third row : [0 A B C 0 ... 0] ;
% [...]
% - Last row = [0 ... 0 1].

% We implement mat : 
Mat = spalloc(N+1, N+1, 3 * (N - 1) + 2); % Sparse matrix.
Mat(1, 1) = 1; % First row.

for i=2:N % Intermediate rows.
    Mat(i, [i-1, i, i+1]) = [A B C];
end

Mat(end, end) = 1; % Last row.
%
%% 3. Implicit Euler - Backward in time loop :
% Vector of price (for all x) that we will update backwardly from T to 0 
% (i.e from tM to t0) :
c = max(S0 * exp(x') - K, 0); % @ Maturity : "v(T, x) = (S0 * e^x - K)^+".

% RHS with the known part of the form "(-1/dt) * V_{j+1, i}" :
rhs = zeros(N+1,1);

for j=M:-1:1 % BACKWARD from t_M to t_0.
    % We know c t_j -> we compute cnew t_{j-1}.
    rhs(1) = 0; % First element [BOUNDARY CONDITION].
    rhs(2:end-1) = -c(2:end-1);
    rhs(end) = Smax - K * exp(-r * (T - (j - 1) * dt)); % Last element [BOUNDARY CONDITION].
    
    % Invert the relation to get V (named c here), where all the "V_{j, i}"
    % terms are (the UNKNOWN !) :
    c = Mat \ rhs;
end
%
%% Plot "call price at t=0" VS "S at t=0" (for all S in [Smin, Smax]) :
% S at t = 0 represents the spot price (so, now). S0 (cf "Model
% parameters") is one example of possible spot price.
figure
plot(S0 * exp(x), c); 
title('Call price at t=0 vs (all possible) Spot Price S_0'); 
xlabel('S at t = 0 (spot price)');
ylabel('Call price at t = 0');
%
%% Compute the price at t = 0 :
price = interp1(x, c, 0, 'spline') % Interpolation of c at x = 0, i.e S = S0
% (S0 being the one we the one we initialized in "Model parameters").
[price_ex, ~] = blsprice(S0, K, r, T, sigma) % For comparison : B&S price.





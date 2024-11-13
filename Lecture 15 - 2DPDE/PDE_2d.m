%% LECTURE 15 - 2D PDE - Noé Debrois - 13/11/2024
% This code implements the pricing of a Knock&Out Call Option on a Basket
% of two Underlying Assets which are NOT CORRELATED. See the other file for
% the case with correlation. We use a PDE method here with discretisation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear;
%
%% Parameters :
% Option parameters :
T = 1;
S10 = 2; S20 = 1; % Initial asset prices
% w1 = w2 = 1/2 here : cf the bottom of this file.
K = 0.9;
D1 = 0.8; U1 = 4; D2 = 0.5; U2 = 2; % Up and Down barriers for the 2 assets
r = 0.001; sigma1 = 0.6; sigma2 = 0.4;

% Space grid parameters :
x1min = log(D1 / S10); x1max = log(U1 / S10);
x2min = log(D2 / S20); x2max = log(U2 / S20);
N1 = 50; N2 = 50; % Number of points in each axis.
x1 = linspace(x1min, x1max, N1 + 1); dx1 = x1(2) - x1(1);
x2 = linspace(x2min, x2max, N2 + 1); dx2 = x2(2) - x2(1);
num_points = (N1 + 1) * (N2 + 1); % Total number of nodes in the grid

% Time grid parameters :
Mt = 100; dt = T/Mt;
%
%% Build the grid :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               i+1
%   i-(N2-1)     i    i+(N2+1)
%               i-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundaries :
% We have to be able to know WHEN we are on the boundary of the domain : we
% have 4 edges that we call West, South, North, East.
W = [1:N2 + 1]; % West boundary
S = [1:N2 + 1:num_points]; % South
N = [(N2 + 1):N2 + 1:num_points]; % North
E = [S(end) + (0:N2)]; % East
Boundary = sort(unique([N E S W]));

% Matrix M is built according to the equation and its discretization :
M = spalloc(num_points, num_points, num_points * 5);
for i=1:num_points
    if min(abs(i-Boundary))==0 % i is a boundary index
        M(i,i) = 1;
    else
        %
        % Look at the coefficients in the discretized equation : here they
        % are implemented. cf bottom of page 67 in my notes (PDE.pdf).
        % For example : M(i, i) contains the coeffs that multiply V_j_i.
        %
        % time derivative
        M(i,i) = -1 / dt;
        % no derivative term
        M(i,i) = M(i,i) - r;
        % first order derivative w.r.t x1
        M(i, i + (N2 + 1)) = (r - sigma1^2 / 2) / (2 * dx1);
        M(i, i - (N2 + 1)) = - (r - sigma1^2 / 2) / (2 * dx1);
        % first order derivative w.r.t x2
        M(i, i + 1) = (r - sigma2^2 / 2) / (2 * dx2);
        M(i, i - 1) = - (r - sigma2^2 / 2) / (2 * dx2);
        % second order derivative w.r.t x1
        M(i, i + (N2 + 1)) = M(i, i + (N2 + 1)) + (sigma1^2 / 2) / dx1^2;
        M(i, i) = M(i, i) - 2 * (sigma1^2 / 2) / dx1^2;
        M(i, i - (N2 + 1)) = M(i, i - (N2 + 1)) + (sigma1^2 / 2) / dx1^2;
        % second order derivative w.r.t x2
        M(i, i + 1) = M(i, i + 1) + (sigma2^2 / 2) / dx2^2;
        M(i, i) = M(i, i) - 2 * (sigma2^2 / 2) / dx2^2;
        M(i, i - 1) = M(i, i - 1) + (sigma2^2 / 2) / dx2^2;
    end
end
%
%% Backward in time procedure
[X1, X2] = meshgrid(x1, x2); % This gives us 2 matrices : coordinates of both x1 and x2.
V = max(0.5 * (S10 * exp(X1(:)) + S20 * exp(X2(:))) - K, 0); % (w1 = w2 = 1/2)
V(Boundary) = 0;
for j=Mt:-1:1
    rhs = - V / dt; 
    rhs(Boundary) = 0; % Boundary Condition
    V = M \ rhs;
end
Vmat = reshape(V, size(X1));
%
%% Plot
figure
surf(X1, X2, Vmat) % Get the price surface
price = griddata(S10 * exp(X1), S20 * exp(X2), Vmat, S10, S20) % Get the price
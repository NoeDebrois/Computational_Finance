%% LECTURE 16 - Heston Model - Noé Debrois - 22/11/2024
% This code implements the QE Scheme for Heston model, from the article of
% Andersen : "Efficient Simulation of the Heston Stochastic Volatility 
% Model".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = Heston_QE(X0, r, x, T, Nsim, N)
%--- PARAMETERS -----------------------------------------------------------
epsilon = x(1); % Vol-of-var.
k = x(2); % Speed of mean reversion speed.
rho = x(3); % Correlation between price and variance.
theta = x(4); % Mean.
dt = T / N; % Time grid.
V0 = x(5); % Variance at time 0.
V = zeros(N + 1, Nsim); % Variance. 
X = zeros(N + 1, Nsim); % Price.

% Random changes for volatility are sampled from one of two distributions,
% depending on the ratio Psi = s^2 / m^2, where m & s are mean and variance
% of next volatility value, conditioned to current one.
% - Scheme 1 (exponential approx.) selected when Psi > Psi_cutoff "E" ;
% - Scheme 2 (quadratic approx.) selected when 0 < Psi < Psi_cutoff "Q".

% Choose Psi_cutoff = 1.5 as suggested in article :
Psi_cutoff = 1.5;

%--- DISCRETIZE V ---------------------------------------------------------
V(1,:) = V0;
for i=1:N
    % STEP 1 & 2 :
    % Calculate m, s, cf slide 18 :
    m = theta + (V(i, :) - theta) * exp(-k * dt);
    m2 = m.^2;
    % cf slide 18 :
    s2 = V(i, :) * epsilon^2 * exp(-k * dt) * (1 - exp(-k * dt)) / k +...
        theta * epsilon^2 * (1 - exp(-k * dt))^2 / (2 * k);
    s = sqrt(s2);

    % Compute Psi :
    Psi = (s2)./(m2);
    
    % STEP 3, 4, 5 :
    % Depending on Psi, use exp (E) or quad (Q) scheme to calculate next V:
    index = find(Psi > Psi_cutoff); % Detection for E scheme.

    % Exponential approximation, used for Psi>Psi_cutoff) : CF SLIDE 23.
    %
    % PDF of V(t+dt) is p * delta(0) + (1 - p) * (1 - exp(-beta * x)
    % thus a prob. mass in 0 & an exp tail after that :
    p_exp = (Psi(index)-1)./(Psi(index)+1);	% Prob. mass in 0, Eq 29 ;
    beta_exp = (1-p_exp)./m(index); % Exponent of exp. density tail, Eq 30;

    % Gets x from inverse CDF applied to uniform U (SLIDE 21) :
    U = rand(size(index));
    V(i + 1, index) = (log((1 - p_exp) ./ (1 - U)) ./ beta_exp) .* (U > p_exp);

    index = find(Psi <= Psi_cutoff); % Detection for Q scheme.

    % Quadratic approximation, used for 0 < Psi < Psi_cutoff : CF SLIDE 23.
    %
    % V(t + dt) = a * (b + Zv)^2, Zv ~ N(0, 1). SLIDE 24.
    % From Slide 23 :
    invPsi = 1 ./ Psi(index);
    % Eq 27 : 
    b2_quad = 2 * invPsi - 1 + sqrt(2 * invPsi) .* sqrt(2 * invPsi - 1);
    % Eq 28 :
    a_quad = m(index) ./ (1 + b2_quad);	
    % Eq.23 (slide 24) :
    V(i + 1, index) = a_quad .* (sqrt(b2_quad) + randn(size(index))) .^2 ;
end

%--- DISCRETIZE X ---------------------------------------------------------
X(1, :) = X0;
% Central discretization scheme
gamma1 = 0.5;
gamma2 = 0.5; 

% Discretization (cf slide 27) :
k0 = r * dt - rho * k * theta * dt / epsilon;
k1 = gamma1 * dt * (k * rho / epsilon - 0.5) - rho / epsilon;
k2 = gamma2 * dt * (k * rho / epsilon - 0.5) + rho / epsilon;
k3 = gamma1 * dt * (1 - rho^2);
k4 = gamma2 * dt * (1 - rho^2);

for i=1:N
    % EQ 33 :
    X(i + 1, :) = exp(log(X(i, :)) + k0 + k1 * V(i, :) + k2 * V(i + 1, :)...
        + sqrt(k3 * V(i, :) + k4 * V(i + 1, :)) .* randn(1, Nsim));
end

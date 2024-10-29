%% LECTURE 7 - Carr-Madan Method - Noé Debrois - 27/10/2024
% This code implements the Carr-Madan (CM) algorithm for pricing a plain 
% vanilla call option using the B&S model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Price] = FFT_CM_Call_BS(K_i)
    % Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm.
    % Model: B&S.
    % For explanation on Carr-Madan (CM): see FFT.pdf.

    % [Price] = FFT_CM_Call_BS(K_i) - could be a vector
    % INPUT : K_i = strikes = could be a vector
    % OUPUT : Price = prices

    % Model parameters :
    S0 = 100;   % Spot price
    T = 1;      % Maturity
    r = 0.0367; % rf interest rate
    sig = 0.17801; % volatility

    % Discretization parameters:
    %
    % 1st step: Truncate the integral from -A/2 to +A/2.
    %
    % - Npow: Truncation of the domain (requires a large domain, hence many points N).
    % - N: Number of grid points (large, but FFT is fast).
    % - A: Upper bound for the integral domain.
    Npow = 15;
    N = 2^Npow;
    A = 600; % Truncation of the integral.

    % Define v_j (or j*eta in the article) to compute the integral as a summation:
    %
    % 2nd step: Quadrature formula to compute z_T(k).
    %
    % - v = v_j = j * eta (article) [WARNING: the article multiplies by 2, 
    %   so j * eta goes from 0 to (N-1) * eta].
    eta = A / N;
    v = [0 : eta : A * (N - 1) / N]; % Trapezoïdal rule with nodes j * eta
    % (j ranges from 0 to N-1).
    v(1) = 1e-22; % Small value to avoid division by zero in FFT{z_T}(0) 
    % calculation. Correction : it can't be equal to zero (otherwise NaN).

    % Define lambda to compute the summation via FFT:
    %
    % 3rd step: Compute F_n via FFT (or DFT).
    %
    % - lambda: Step size in log-strike grid.
    lambda = 2 * pi / (N * eta);
    k = -lambda * N / 2 + lambda * (0:N-1); % log strike grid (cf article).
    %
    %% Pricing :
    tic

    % Fourier Transform of z_k (Z_k = F(z_k)) :
    Z_k = trans_fourier_BS(r, sig, T, v);
    
    % Trapezoïdal quadrature formula :
    w = ones(1, N);
    w(1) = 0.5;
    w(end) = 0.5;

    % Argument inside FFT():
    x = w .* eta .* Z_k .* exp(1i * pi * (0:N-1));

    % FFT formula for z_T(k_l): 
    %
    % ! fft() because MATLAB's FFT is implemented with a "-" in the exp,
    % whereas probabilists use the IFT with "-". So here we perform the FFT
    % in MATLAB's POV, but the IFT in probabilists' POV.
    %
    z_k = real(fft(x) / pi); % Due to numerical approximation, use 'real()'
    % to put the imaginary part to zero.
    % We know theoretically that it should be indeed real.

    % Compute option prices ; C(k) = z_T(k) + (1 - exp(k - r * T))^+.
    C = S0 * (z_k + max(1 - exp(k - r * T), 0)); % Option prices array.

    time = toc % print the time for pricing.
    %
    %% Post-processing :

    % - Convert log-strikes to strikes:
    K = S0 * exp(k);

    % Output Processing:
    % - Filter strikes: remove very small and large strikes.
    index = find(K > 0.1 * S0 & K < 3 * S0);
    C = C(index);
    K = K(index);
    %
    %% Plot the results :
    % - WARNING: This is not a put option, as the x-axis represents the 
    % strike, not the spot price.
    figure
    plot(K, C, 'r');
    hold on
    axis([0 2*S0 0 S0]);
    title('Option Price vs Strike');
    ylabel('Option Price');
    xlabel('Strike');
    %
    %% Interpolation :
    Price = interp1(K, C, K_i, 'spline');
end

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function zeta_T = trans_fourier_BS(r, sig, T, v)
    % Computation of Zeta_T(v) [article] = g_T(v) [lecture notes].
    % Cf formula (11.19) in the lecture notes, page 15 :
    % zeta_T(v) = exp(ivrT) * [Phi_T(v-i) - 1] / [iv * (1+iv)]
    %
    zeta_T = (exp(1i * v * r * T)) .* ...
        ((characteristic_func_BS(sig, T, v - 1i) - 1) ./ ...
        (1i * v .* (1 + 1i * v)));
    
end

function phi_T = characteristic_func_BS(sig, T, y)
    % Risk-Neutral characteristic function computation.
    % Formula for a Geometric Brownian Motion (GBM) :
    % phi_T(y) = exp(- sigma^2/2 * iyT - sigma^2/2 * y^2 * T)
    %
    phi_T = exp(1i * (-sig^(2) / 2 * T) .* y ...
        - (T * sig.^(2) * y.^(2)) / 2);
end

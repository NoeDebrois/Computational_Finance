%% LECTURE 10 - Convolution / FFT Methods - Noé Debrois - 30/10/2024
% This code implements the algorithm that uses "kernel.m" and goes backward
% in time, according to the algorithm seen in class (and the one from the
% article CONV_Method.m).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,v] = CONV(S_0, K, Ndate, N, Barrier, param)
%
b = 2.5; % xmax
[x, ~, ~, H] = kernel(N, -b, b, param, 0);
S = S_0 * exp(x); % Spot price on the log-price grid.

% Payoff and Fourier transform of the density :
v = max(S - K, 0) .* (S > Barrier); % Down & Out here ! 
H = ifftshift(H); % Once again we have to transform H to be in Matlab style

for j = Ndate:-1:1 % Backward in time procedure :
    % From T to 0 (T, T-delta, T-2delta, ..., 0) :
    % 1) Take the FT of the payoff
    % 2) Take the FT of the density = characteristic function (conjugated!)
    % 3) Multiply the two and multiply by exp(-r*delta) (discounting)
    % 4) Take the IFT
    % 5) Put everything under the barrier to 0 (of above if Up&Out)
    % And repeat the procedure till t = 0.
    %
    % C(T-delta) = IFT(exp(-r * delta) * FT(payoff) * Phi_T^*)
    %
    v = real((fft(ifft(v) .* H))) * exp(-param.rf * param.dt);
%   Another way to write it :    
%   v = real(fftshift(fft ...
%         (ifft(ifftshift(v)) .* H))) * ...
%         exp(-param.rf * param.dt);
    v(S <= Barrier) = 0; % Down & Out here ! If we remove this, we would 
    % a classical EU option.
end

index = find((S > 0.1 * S_0) .* (S < 3 * S_0));
S = S(index); v = v(index);
figure
plot(S, v);
xlabel('Grid of spot price S'); 
ylabel('Call Option Price (v)')
title('Barrier (Down&Out) Call Option Price (v) vs a grid of Spot Price (S)');
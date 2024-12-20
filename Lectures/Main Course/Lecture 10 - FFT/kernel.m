%% LECTURE 10 - Convolution / FFT Methods - No� Debrois - 30/10/2024
% This code implements the computations of the Complex Conjugate of the
% Characteristic Function of X, on the Fourier space grid, and the Density
% Function of -X on the Log-price grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, h, w, H] = kernel(ngrid, xmin, xmax, parameters, flag)
% cf FFT.pdf & CONV.pdf : in the computations, to make a convolution appear
% we need to complex-conjugate the characteristic function. That's why :
% flag = 0 for complex conjugate -> characteristic fct backward-in-time ; 
% flag = 1 for NO complex conjugate -> characteristic fct forward-in-time.

if nargin == 4
    flag = 0; % By default, we need the complex conjugate.
end
    
% We need two grids :
%% Log-price grid :
N = ngrid / 2;
dx = (xmax - xmin) / ngrid;
x = dx * (-N:N-1);
%
%% Fourier space grid :
dw = 2 * pi / (xmax - xmin);
w = dw * (-N:N-1);
%
%% Compute conjugated char fct (H) on Fourier grid, & density (h) of -X :
H = charfunction(w, parameters, flag); %  H = Phi_X^*, conjugated char fct 
% of X (cf charfunction.m file).
%
% We need to correct the TWO shifts (both grids) induced by Matlab :
h = real(fftshift(fft(ifftshift(H)))) / (xmax - xmin); % h is the density 
% of (-X) (since H = Phi_X^*). "fft" is the MATLAB's fft, i.e, the
% probabilists' IFT (for probabilists, PHI = FT(density)).
%
% Explanation of the shifts for h, from left to right :
% - fftshift() : this corrects the log-price grid (x) which doesn't go from
% 0 to 2N BUT from -N to N-1 ;
% - ifftshift() : this corrects the Fourier grid (w) which doesn't go from
% 0 to 2N BUT from -N to N-1.
%
% cf comments at the end of this file.
%
%% Plots :
% Density of -X :
figure
plot(x, h)
title("Density of -X, on the log-price grid (x)");
xlabel("Log-price grid (x)");
ylabel("Density of -X")

% Conjugated characteristic function :
figure
plot(w, H)
title("Conjugated characteristic function of X, on the Fourier space grid (w)");
xlabel("Fourier space grid (w)");
ylabel("Conjugated characteristic function of X")
%
%% Explanations :
% Basically, MATLAB implements an fft (as do other sets of functions) 
% whose results AND inputs are ordered (for an array of N elements), 
% from 0 to (N/2-1) and then from �N/2 to -1. 
% (I will call this 'swapped' from now on.) 
% But in what I might call a 'natural' ordering, 
% for a spectrum or time function centred at zero and extending equal time 
% (or frequency) either side of zero, 
% the arrays of input data and results would be ordered (2�s complement style)
% from �N/2 to N/2-1 where N is the number of elements. 
% That is, time zero (or frequency zero) are at the centre of such an array. 
% 
% It is well known that the results of MATLAB's fft() function must be 
% re-ordered with the MATLAB function fftshift() so that they are 'naturally' 
% ordered (for instance when plotting).
% 
%    H = fftshift( fft ( h ) );
% 
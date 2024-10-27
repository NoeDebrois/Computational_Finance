%% LECTURE 7 - Carr-Madan Method - Noé Debrois - 27/10/2024
% This code is calling the Carr-Madan (CM) algorithm for pricing a plain 
% vanilla call option using the Kou model for jump-diffusion processes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
% Parameters :
Strike = [80 90 100 110];
S0 = 102;
T = 1;
r = 0.01 / 100;
params = [0.5 3 0.6 20 30];
% Call : 
FFT_CM_Call_Kou(Strike, params, T, r, S0)
%
% It's very fast and with a nice precision.
% It should be great to calibrate to the market, because calibration is a 
% process that requires a lot of iterations, but on plain vanilla objects.
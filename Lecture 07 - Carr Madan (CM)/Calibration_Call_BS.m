%% LECTURE 7 - Carr-Madan Method - Noé Debrois - 27/10/2024
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
%% Data pre-processing :

% European Call Options : we use data with two maturities, and a few
% strikes for each maturity.

Maturity = 4/12; % Maturity of 4 months 
Data = [25.30	200 % Prices, Strikes
        10.00	225
        18.20	210
        41.85	180
        15.20	215
        7.85	230
        12.36	220
        46.05	175
        70.70	150
        50.85	170
        60.40	160
        29.40	195];
Data1 = [Data(:,2), Maturity * ones(size(Data(:,2))), Data(:,1)];
% Data1 is in the form : [Strike, Time to Maturity, Market Price].

Maturity = 4/252; % Maturity of 4 (trading) days
Data = [0.24	230
        9.25	210
        2.26	220
        0.76	225
        19.03	200
        5.20	215
        14.00	205
        0.10	235];
Data = [Data(:,2), Maturity * ones(size(Data(:,2))), Data(:,1)];
% Data is in the form : [Strike, Time to Maturity, Market Price].

Data = [Data1; Data]; % We merge the two datasets (with two different 
% maturities).
%
%% Optimisation :

% WE WANT to compute the sigma that minimizes the distance between market
% price and model price, where the model price is computed using fun.m (i.e
% Kou model, by Carr-Madan algorithm).

% Parameters :
r = 0.02; spot = 218.75;

% Optimisation :
x0 = [0.5, 5, 0.8, 15, 15]; % Initial guess
LB = [0.1, 0, 0, 1,1];      % Lower bound
UB = [0.8, 20, 1, 50, 50];  % Upper bound

tic

[params, error] = lsqnonlin(@(params) ...
        fun_CM_BS(params, spot, Data(:,1), r, Data(:,2), Data(:,3)), ...
        x0, LB, UB)

toc
%
%% Plot :

figure
plot(Data(:,1), Data(:,3), '+', 'markersize', 16, 'linewidth', 5); 
hold on
% Compute and plot the calibrated prices (so with the calibrated sigma) :
Price = fun(params, spot, Data(:,1), r, Data(:,2), 0);
plot(Data(:,1), Price, 's', 'markersize', 16, 'linewidth', 5); 
legend('Market price', 'Model price')
xlabel("Strike")
ylabel("Price")
title("Calibration using CM algorithm, B&S model")
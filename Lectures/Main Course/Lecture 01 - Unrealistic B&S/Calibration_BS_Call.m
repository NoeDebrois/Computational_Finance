% Calibration of the volatility to find the best sigma, i.e, the sigma such
% that the B&S price is as close as possible to the market price.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
% European Call Options
Maturity = 4/12; % 4 months
% Format of the data : Price, Strike
Data = [25.30	200  
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
% Format of the data : Strike, Maturity, Price
Data1 = [Data(:,2), Maturity * ones(size(Data(:,2))), Data(:,1)];
Data = Data1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters :
r = 0.02; % Risk-free rate
spot = 218.75; % Spot price (S0)
%
% Calibrate the volatility, 'vol', such that the B&S price is as close as
% possible to the market price (in variable Data). 
% -> 'fun' (cf fun.m) computes the difference between the B&S price (B&S
% formula) and the real market price (in variable Data) :
[sigma, error] = lsqnonlin(@(vol)fun(vol, spot, Data(:,1), r, Data(:,2), Data(:,3)), ...
    0.2, 0.01, 0.8) % x_0, lower_bound, upper_bound (for the volatility).
% The variable 'sigma' contains the calibrated volatility value.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the result of the calibration :
figure
plot(Data(:,1),Data(:,3),'+'); hold on
% B&S price with the calibrated volatility 'sigma' :
Price = fun(sigma, spot, Data(:,1), r, Data(:,2), 0); 
% Market price :
plot(Data(:,1), Price, 's'); 
legend('Market price', 'BS price')
%  Here we show that B&S model is unrealistic by plotting the volatility
%  smile obtained through B&S model. We clearly see that it's not constant
%  contrary to what we assume in B&S model :
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
Data = [Data(:,2), Maturity * ones(size(Data(:,2))), Data(:,1)];
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters :
r = 0.02; % Risk-free interest rate
spot = 218.75; % Sport price (S0)
sigma = 0.2; % volatility
% 
options = optimoptions('fsolve', 'FunctionTolerance', 1e-10, 'StepTolerance', 1e-8, 'OptimalityTolerance', 1e-8, 'MaxIterations', 800, 'MaxFunctionEvaluations', 500)
% Invert the B&S formula to find the volatility that permits to have
% equality between B&S price and Market Price :
for i=1:size(Data,1)
    sigma = fsolve(@(vol)fun(vol, spot, Data(i,1), r, Data(i,2), Data(i,3)), sigma, options);
    Vol(i) = sigma;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the volatility smile we obtain :
figure
plot(Data(:,1), Vol, 'd'); 
xlabel('Strike'); 
ylabel('Implied volatility');
title("Volatility smile obtained by inversion of the B&S formula")
% The volatility looks like a smile : it is not constant.
% This shows that we need better models than the B&S model !
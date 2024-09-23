% Here we show that B&S model is unrealistic by plotting the market log-
% returns of Apple vs normally-distributed log-returns. We also show that
% the skewness and kurtosis of the market log-returns are not consistent
% with the skewness and kurtosis of a normal distribution :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;
%
Stock % Upload Data from Apple stock.
Apple_MarketData = flipud(Apple_MarketData); % First element is the oldest.
% Compute the (daily-)log returns (cf Ptf management course) :
logreturn = log(Apple_MarketData(2:end) ./ Apple_MarketData(1:end-1));
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot of the (daily-)log returns from Apple :
figure; hold on;
plot(logreturn); % Plot the time-series of (daily-)log returns.
title("Apple's (daily-)log returns time series");
figure;
hist(logreturn, 40); % 2nd argument : number of bins.
title("Apple's (daily-)log returns histogram");
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of B&S parameters from a historical time series :
mean_value = mean(logreturn);
variance = var(logreturn);
dt = 1/252; % 252 trading days per year.
% cf lecture 1 of Computational Finance (intro) :
sigma = sqrt(variance / dt)
mu = mean_value / dt + sigma^2 / 2
% Check Normality assumption :
figure
qqplot(logreturn)
sk = skewness(logreturn) % The skewness for a normal distribution is zero, 
% and any symmetric data should have a skewness near zero. 
ku = kurtosis(logreturn) % Normal distribution has a kurtosis = 3.
% The market log-returns are not completely normally-distributed...
% This shows that we need better models than the B&S model !
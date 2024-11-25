%% Timing strategy with moving averages - Lesson 1.

clear all
close all
clc
%% Read Prices

path_map = '/Users/noedebrois/Desktop/Desktop - Noé’s MacBook Air/Politecnico/Computational Finance/AY 2024:2025/Portfolio Management/';
filename = 'spx_price.xlsx';
% Read the excel table :
table_prices = readtable(strcat(path_map, filename));
%% Transform prices from table to timetable

dt = table_prices(:,1).Variables;
values = table_prices(:,2).Variables;
nm = table_prices.Properties.VariableNames(2);

% "timetable" is a type of table that associates a time with each row for 
% use with time series data. So we convert the data to a "timetable" :
myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', nm); 
%% Building the Strategy

closePrices = myPrice_dt.Variables;
dates = myPrice_dt.Time;

% Defining window size
windowShort = 50;  % short
windowLong = 200;  % long

% Compute moving average
SMA50 = movmean(closePrices, windowShort);
SMA200 = movmean(closePrices, windowLong);

% Initialize signal (buy/sell signal)
signals_ = zeros(length(closePrices), 1);

% Create signals :
% - A 1 (buy signal) is generated if the short moving average (SMA50) 
% crosses above the long moving average (SMA200).
% - A 0 (sell/exit signal) is generated if the short moving average crosses
% below the long moving average.
for i = windowLong:length(closePrices)
    if SMA50(i) > SMA200(i)
        signals_(i) = 1;  % buy signal
    elseif SMA50(i) < SMA200(i)
        signals_(i) = 0; % out signal
    else
        signals_(i) = 0;  % no action
    end
end

% Plot 1
figure;
% Plot S&P500 index :
plot(dates, closePrices, 'k', 'DisplayName', 'Closing price');
hold on;
% Plot MA on a 50 days window :
plot(dates, SMA50, 'b', 'DisplayName', 'SMA 50 days');
hold on;
% Plot MA on a 200 days window :
plot(dates, SMA200, 'r', 'DisplayName', 'SMA 200 days');
legend('SPX Index', 'MA 50', 'MA 200')
% Add Buy/Sell signals on plot :
buySignals = signals_ == 1;
sellSignals = signals_ == 0;
% Add green upward triangles for buy signals (g^) :
plot(dates(buySignals), closePrices(buySignals), 'g^', 'MarkerSize', 8, 'DisplayName', 'BUY');
% Add red downward triangles for sell signals (rv) :
plot(dates(sellSignals), closePrices(sellSignals), 'rv', 'MarkerSize', 8, 'DisplayName', 'SELL');
title('Moving Average Crossover Strategy on S&P 500');
xlabel('Date');
ylabel('Closing Price');
legend('show');
grid on;
hold off;
%% Equity

% Compute Equity Curve
index = windowLong;
% Daily returns from the closing prices :
ret = closePrices(index:end,:)./closePrices(index-1:end-1,:);
Signal = signals_(index:end,:);
% Returns based on the signals generated :
ret_strategy = Signal(1:end-1).*ret(2:end,:); % return of the strategy

InitialCapital = 100;
% Initialise the equity curve (+1 to include the initial capital) :
EquityCurve = zeros(length(ret_strategy)+1, 1); 
EquityCurve(1) = InitialCapital;

% For loop to compute the Equity Curve (cf Lecture 1).
% It tracks the portfolio's growth over time, adjusting only when a buy
% signal is active :
for t = 2:length(EquityCurve)
    if Signal(t-1) == 1
        % If the signal is 1, apply the return :
        EquityCurve(t) = EquityCurve(t-1) * (ret_strategy(t-1));
    else
        % If the signal is 0, keep the equity unchanged :
        EquityCurve(t) = EquityCurve(t-1);
    end
end

closePrices_subset = closePrices(index:end);
closePrices_subset = 100*closePrices_subset/closePrices_subset(1);

% Plot 2
% The cumulative performance of the trading strategy (EquityCurve) compared
% to the normalized S&P 500 price in the same period.
figure()
title("Cumulative performance of the trading strategy VS normalized S&P500 price");
yyaxis left;
plot(dates(index:end), EquityCurve, 'r', 'linewidth', 4)
hold on
plot(dates(index:end), closePrices_subset, 'k')
yyaxis right;
plot(dates(index:end), Signal(1:end), 'b')
legend('Equity Curve', 'Normalized S&P500 price', 'BUY/SELL signal')

% Last 2 years
EquityCurve_sub = EquityCurve(length(EquityCurve)-500:end);
EquityCurve_sub = 100*EquityCurve_sub/EquityCurve_sub(1);
closePrices_subset = closePrices_subset(length(closePrices_subset)-500:end);
closePrices_subset = 100*closePrices_subset/closePrices_subset(1);

% Plot 3
% Same as plot 3 but only on the last 2 years.
figure()
plot(dates(length(dates)-500:end), EquityCurve_sub, 'r', 'linewidth', 4)
hold on
plot(dates(length(dates)-500:end), closePrices_subset, 'k')
legend('Equity Curve', 'Normalized S&P500 price')
title("Performance of the strategy in the last two years");

% The equity curve doesn't grow as fast as the S&P price curve. So the
% strategy is underperforming compared to the S&P index. It means that the
% strategy generates less returns than simply holding the S&P 500 index 
% over the same period. But using the crossover in MAs we try to avoid big
% losses so it's probably quite great to reduce risk.
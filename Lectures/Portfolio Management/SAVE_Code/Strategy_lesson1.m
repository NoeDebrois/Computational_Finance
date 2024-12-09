clear all
close all
clc

%% Read Prices
path_map        = folder_path;
filename        = 'spx_price.xlsx';

table_prices = readtable(strcat(path_map, filename));
%% Transform prices from table to timetable

dt = table_prices(:,1).Variables;
values = table_prices(:,2).Variables;
nm = table_prices.Properties.VariableNames(2);

myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', nm); %già in datetime
%% Building the Strategy
closePrices = myPrice_dt.Variables;
dates = myPrice_dt.Time;

% Defining window size
windowShort = 50;  % short
windowLong = 200;  % long

% Compute mov average
SMA50 = movmean(closePrices, windowShort);
SMA200 = movmean(closePrices, windowLong);

% Initialize signal
signals_ = zeros(length(closePrices), 1);

% Create signals
for i = windowLong:length(closePrices)
    if SMA50(i) > SMA200(i)
        signals_(i) = 1;  % buy signal
    elseif SMA50(i) < SMA200(i)
        signals_(i) = 0; % out signal
    else
        signals_(i) = 0;  % No action
    end
end

% Plot
figure;
plot(dates, closePrices, 'k', 'DisplayName', 'Prezzo di chiusura');
hold on;
plot(dates, SMA50, 'b', 'DisplayName', 'SMA 50 giorni');
hold on;
plot(dates, SMA200, 'r', 'DisplayName', 'SMA 200 giorni');
legend('SPX Index', 'MA 50', 'MA 200')
% Add Buy/Sell signals on plot
buySignals = signals_ == 1;
sellSignals = signals_ == 0;

plot(dates(buySignals), closePrices(buySignals), 'g^', 'MarkerSize', 8, 'DisplayName', 'Acquisto');
plot(dates(sellSignals), closePrices(sellSignals), 'rv', 'MarkerSize', 8, 'DisplayName', 'Vendita');
title('Moving Average Strategy on S&P 500');
xlabel('Date');
ylabel('Closing Price');
legend('show');
grid on;
hold off;
%% Equity
% Compute Equity Curve
index = windowLong;
ret = closePrices(index:end,:)./closePrices(index-1:end-1,:);
Signal = signals_(index:end,:);
ret_strategy = Signal(1:end-1).*ret(2:end,:);

InitialCapital = 100;
EquityCurve = zeros(length(ret_strategy)+1, 1);  % +1 per includere il capitale iniziale
EquityCurve(1) = InitialCapital;

% Ciclo per calcolare la curva di equity
for t = 2:length(EquityCurve)
    if Signal(t-1) == 1
        % Se il segnale è 1, applica il rendimento
        EquityCurve(t) = EquityCurve(t-1) * (ret_strategy(t-1));
    else
        % Se il segnale è 0, mantieni l'equity invariata
        EquityCurve(t) = EquityCurve(t-1);
    end
end

closePrices_subset = closePrices(index:end);
closePrices_subset = 100*closePrices_subset/closePrices_subset(1);

% Plot
figure()
yyaxis left;
plot(dates(index:end), EquityCurve, 'r', 'linewidth', 4)
hold on
plot(dates(index:end), closePrices_subset, 'k')
yyaxis right;
plot(dates(index:end), Signal(1:end), 'b')

mean(Signal)
% Last 2 years
EquityCurve_sub = EquityCurve(length(EquityCurve)-500:end);
EquityCurve_sub = 100*EquityCurve_sub/EquityCurve_sub(1);
closePrices_subset = closePrices_subset(length(closePrices_subset)-500:end);
closePrices_subset = 100*closePrices_subset/closePrices_subset(1);

% Plot
figure()
title('Performance of Strategy in the last two years')
plot(dates(length(dates)-500:end), EquityCurve_sub, 'r', 'linewidth', 4)
hold on
plot(dates(length(dates)-500:end), closePrices_subset, 'k')

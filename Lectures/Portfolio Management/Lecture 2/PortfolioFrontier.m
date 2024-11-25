clear all
close all
clc

%% Read Prices
path_map        = '/Users/noedebrois/Desktop/Desktop - Noé’s MacBook Air/Politecnico/Computational Finance/Computational_Finance/Portfolio Management/Lecture 2/';
filename        = 'asset_prices_lesson2.xlsx';

table_prices = readtable(strcat(path_map, filename));
%% Transform prices from table to timetable
dt = table_prices(:,1).Variables;
values = table_prices(:,2:end).Variables;
nm = table_prices.Properties.VariableNames(2:end);
myPrice_dt = array2timetable(values, 'RowTimes', dt,'Variablenames', nm); 
%% Selection of a subset of Dates
start_dt = datetime('01/02/2021', 'InputFormat', 'dd/MM/yyyy');
end_dt = datetime('01/04/2021', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt, 'closed');
subsample = myPrice_dt(rng,:);
prices_val = subsample.Variables;
dates_ = subsample.Time;
%% Calculate returns
% Method 1 
ret = prices_val(2:end,:)./prices_val(1:end-1,:);
LogRet = log(ret);

% method 2 
LogRet1 = tick2ret(prices_val, 'Method', 'Continuous');
%% Calculate Variance
ExpRet = mean(LogRet);
var_ = var(LogRet);
std_ = std(LogRet);

V = cov(LogRet);
%% Creation of N random portfolio
N = 100000;
NumAssets = 15;
RetPtfs = zeros(1, N);
VolaPtfs = zeros(1, N);
SharpePtfs = zeros(1,N);

for n = 1:N
    w = rand(1, NumAssets);
    w_norm = w./sum(w);
    
    exp_ret_ptf = w_norm*ExpRet';
    exp_vola_ptf = sqrt(w_norm*V*w_norm');
    sharpe_ratio = exp_ret_ptf/exp_vola_ptf;
    
    RetPtfs(n) = exp_ret_ptf;
    VolaPtfs(n) = exp_vola_ptf;
    SharpePtfs(n) = sharpe_ratio;
end
%% Plot
h = figure();
title('Expected Return vs Volatility')
scatter(VolaPtfs, RetPtfs, [], SharpePtfs, 'filled')
colorbar
xlabel('Volatilty')
ylabel('Expected Return')
%% Portfolio Frontier
fun = @(x)x'*V*x;
ret_ = linspace(min(RetPtfs), max(RetPtfs), 100);
x0 = rand(1, NumAssets)';
x0 = x0/sum(x0);

lb = zeros(1, NumAssets);
ub = ones(1, NumAssets);

FrontierVola = zeros(1, length(ret_));
FrontierRet = zeros(1, length(ret_));

for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1, NumAssets); ExpRet];
    beq = [1; r];
    w_opt = fmincon(fun, x0, [], [], Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt'*V*w_opt);
    FrontierVola(i) = min_vola;
    FrontierRet(i) = r;
end
%% Plot
h = figure;
title('Expected return vs volatility')
scatter(VolaPtfs, RetPtfs)
hold on
plot(FrontierVola, FrontierRet)
xlabel('Volatility')
ylabel('Expected return')


%% Add Costraints: min exp = 0.01 and max exp = 0.7 for each asset 

%% Add Costraints - Part 2 : AAPL, AMZN and GOOGL exposition has to be more than 10%


%% Plot the three frontiers


%% Portfolio Frontier respect to benchmark


% Plot


%% Portfolio object

%% COMPUTE EFFICIENT FRONTIER AND PLOT INFORMATION RATIO (IRR)

%% Add turnover costraints and Transaction costs
BuyCost = 0.0020;
SellCost = 0.0020;
Turnover = 0.3; %average turnover < 30%


% Plot

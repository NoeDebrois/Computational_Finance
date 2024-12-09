
clear all
close all
clc

%% Read Prices & Names
load array_prices.mat
load myPrice_dt.mat
start_dt = datetime('01/01/2021', 'InputFormat', 'dd/MM/yyyy'); % dt_4(1)+20;
end_dt   = datetime('01/06/2022', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed'); %fai vedere openLeft ecc 'openLeft'
subsample = myPrice_dt(rng,:);

prices_val = subsample.Variables;
dates_ = subsample.Time;
%% Compute Moments
LogRet = tick2ret(prices_val, 'Method', 'Continuous');
ExpRet = mean(LogRet);
CovMatrix = cov(LogRet);
%% VaR - Historical Simulation - h = 1
% We compute VaR & ES for Equally Weighted Portfolio
weights_EW = ones(size(prices_val, 2), 1).*1/size(prices_val,2);
pRet = weights_EW'*LogRet';
ConfLevel = [0.95, 0.99];

VaR_95 = quantile(pRet, 1-ConfLevel(1,1));
VaR_99 = quantile(pRet, 1-ConfLevel(1,2));
% Expected Shortfall
ES_95 = mean(pRet(pRet < VaR_95));
ES_99 = mean(pRet(pRet < VaR_99));
% Plot
f = figure();
histogram(pRet);
hold on
xline(VaR_95, 'LineWidth', 4, 'Color', 'r')
hold on
xline(VaR_99, 'LineWidth', 4, 'Color', 'm')
legend('Profit & Loss Distr', 'VaR 95%', 'VaR 99%')

%% VaR - Normal distribution: h = 1
mu = mean(pRet);
std_ = std(pRet);

VaR_95Norm = mu-std_*norminv(0.95);
VaR_99Norm = mu-std_*norminv(0.99);
% Matlab Function portvrisk
ValueAtRisk95 = portvrisk(mu, std_, 1-ConfLevel(1,1));
ValueAtRisk99 = portvrisk(mu, std_, 1-ConfLevel(1,2));
% Expected Shortfall
ES_95 = mean(pRet(pRet < VaR_95Norm));
ES_99 = mean(pRet(pRet < VaR_99Norm));

% Plot
histogram(pRet);
hold on
xline(VaR_95Norm, 'LineWidth', 4, 'Color', 'r')
hold on
xline(VaR_99Norm, 'LineWidth', 4, 'Color', 'm')
legend('Profit & Loss Distr', 'VaR 95%', 'VaR 99%')

%% Optmization con target value VaR
fun = @(x) -(ExpRet*x);
x0 = rand(size(LogRet,2),1);
x0 = x0./sum(x0);
lb = zeros(1, size(LogRet,2))+0.001;
ub = ones(1, size(LogRet,2));
Aeq = ones(1, size(LogRet,2));
beq = 1;
pval = 0.05;
tgtVaR = -0.015;

wVar_95 = fmincon(fun, x0, [],[], Aeq, beq, lb, ub, @(x) nonlinConstrVar(x, LogRet, pval, tgtVaR));

VaR_95_test = quantile(wVar_95'*LogRet',0.05);
%% Optimization with ES
tgtES = -0.05;
wES_95 = fmincon(fun, x0, [],[], Aeq, beq, lb, ub, @(x) nonlinConstrES(x, LogRet, pval, tgtES));

VaR_95_TestES = quantile(wES_95'*LogRet', 0.05);
rfPf95 = wES_95'*LogRet';
ES_test = mean(rfPf95(rfPf95<VaR_95_TestES ));

%% Equity Comparison
ret = prices_val(2:end,:)./prices_val(1:end-1,:);
% VaR 95
equityVaR95 = cumprod(ret*wVar_95);
equityVaR95 = 100.*equityVaR95/equityVaR95(1);
[annRet_VaR95, annVol_VaR95, Sharpe_VaR95, MaxDD_VaR95, Calmar_VaR95] = getPerformanceMetrics(equityVaR95);
perfTable_VaR95 = table(annRet_VaR95, annVol_VaR95, Sharpe_VaR95, MaxDD_VaR95, Calmar_VaR95, 'VariableNames',["AnnRet", "AnnVol", "Sharpe", "MaxDD","Calmar"]);


% ES 95
equityES95 = cumprod(ret*wES_95);
equityES95 = 100.*equityES95/equityES95(1);
[annRet_ES95, annVol_ES95, Sharpe_ES95, MaxDD_ES95, Calmar_ES95] = getPerformanceMetrics(equityES95);
perfTable_ES95 = table(annRet_ES95, annVol_ES95, Sharpe_ES95, MaxDD_ES95, Calmar_ES95, 'VariableNames',["AnnRet", "AnnVol", "Sharpe", "MaxDD","Calmar"]);

% EW
wEW = 1/15*ones(15, 1);
equity_ew = cumprod(ret*wEW);
equity_ew = 100.*equity_ew/equity_ew(1);
[annRet_ew, annVol_ew, Sharpe_ew, MaxDD_ew, Calmar_ew] = getPerformanceMetrics(equity_ew);
perfTable_ew = table(annRet_ew, annVol_ew, Sharpe_ew, MaxDD_ew, Calmar_ew, 'VariableNames',["AnnRet", "AnnVol", "Sharpe", "MaxDD","Calmar"]);

% Plot
f1 = figure();
plot(dates_(2:end,1),equityVaR95, 'LineWidth', 2)
hold on
plot(dates_(2:end,1),equityES95, 'LineWidth', 2)
hold on
plot(dates_(2:end,1),equity_ew, 'LineWidth', 2)
legend('Equity VaR', 'Equity ES', 'Equity EW')
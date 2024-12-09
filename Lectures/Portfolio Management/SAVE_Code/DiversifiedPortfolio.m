clear all
close all
clc

%% Read Prices & Names
load array_prices.mat
load myPrice_dt.mat
start_dt = datetime('01/01/2019', 'InputFormat', 'dd/MM/yyyy'); % dt_4(1)+20;
end_dt   = datetime('01/06/2019', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed');
subsample = myPrice_dt(rng,:);

prices_val = subsample.Variables;
dates_ = subsample.Time;
%% Compute Moments
LogRet = tick2ret(prices_val, 'Method','Continuous');
ExpRet = mean(LogRet);
CovMatrix = cov(LogRet);

%% Equally Weighted Ptf
wEW = 1/size(LogRet,2)*ones(15, 1);
RetPtfEW = wEW'*ExpRet';
VolaPtfEW = sqrt(wEW'*CovMatrix*wEW);

%% Risk Contributions of EW Ptf
[relRC_ew, RC_ew, mVol_ew] = getRiskContributions(wEW, LogRet);
% Plot Risk Contribution
pie(relRC_ew, assetNames)

%% Risk Parity
Target = wEW;
Aeq = ones(1,15);
beq = 1;
lb = zeros(1,15);
ub = ones(1,15);
x0 = wEW;
w_RP = fmincon(@(x) mse_risk_contribution(x, LogRet, Target), x0, [], [], Aeq, beq, lb, ub);

[relRC_rp, RC_rp, mVol_rp] = getRiskContributions(w_RP, LogRet);

% Plot
h = figure();
title('Relative Risk Contributions for Risk Parity')
pie(relRC_rp, assetNames)

g = figure();
title('Optimized Weights in Risk Parity Portfolio')
pie(w_RP, assetNames)

%% Most Diversified Portfolio
DR_ew = getDiversificationRatio(wEW, LogRet);
DR_rp = getDiversificationRatio(w_RP, LogRet);

% Maximization of DR
Aeq = ones(1,15);
beq = 1;
lb = zeros(1,15);
ub = ones(1,15);
x0 = wEW;
[w_DR, fval] = fmincon(@(x) -getDiversificationRatio(x, LogRet), x0, [], [], Aeq, beq, lb, ub);

MaxDR = -fval;
[relRC_DR, RC_DR, mVol_DR] = getRiskContributions(w_DR, LogRet);

% Plot Risk Contributions
h = figure();
title('Relative Risk Contributions for Max Div Ratio')
pie(relRC_DR, assetNames)
% Plot Weights
g = figure();
title('Relative Risk Contributions for Max Div Ratio')
pie(w_DR, assetNames)

%% Max Entropy (weights)
EntropyEW = getEntropy(wEW);
EntropyRP = getEntropy(w_RP);
EntropyDR = getEntropy(w_DR);

% Optimization
Aeq = ones(1,15);
beq = 1;
lb = zeros(1,15);
ub = ones(1,15);
x0 = zeros(15,1);
x0(1,1) = 1;
w_MaxEntropy = fmincon(@(x) -getEntropy(x), x0, [], [], Aeq, beq, lb, ub); %from theory
entropy_temp = getEntropy(w_MaxEntropy);

% Max Entropy in Risk Contributions
w_MaxEntropyRisk = fmincon(@(x) -getEntropy(getRiskContributions(x,LogRet)), x0, [], [], Aeq, beq, lb, ub); %from theory
MaxEntropy_RC = getEntropy(w_MaxEntropyRisk);

% Max Entropy in Volatilities
w_MaxEntropyVol = fmincon(@(x) -getEntropy(getVolContributions(x, LogRet)), x0, [], [], Aeq, beq, lb, ub); %from theory
MaxEntropy_Vol = getEntropy(w_MaxEntropyVol);

%% Minimum Variance Portfolio
p = Portfolio('AssetList',assetNames);
p = setDefaultConstraints(p);
PortMV = estimateAssetMoments(p, LogRet,'missingdata',false);
% COMPUTE EFFICIENT FRONTIER
weightsMVP = estimateFrontier(PortMV, 100);
[pfMVP_Risk, pfMVP_Retn] = estimatePortMoments(PortMV, weightsMVP); 

[min_vola, indexMV] = min(pfMVP_Risk);
weightsMV = weightsMVP(:, indexMV);
%% Compare Portfolios (EW, RP, MD, EntropyVol) - In Sample
% Equity Curve
ret = prices_val(2:end,:)./prices_val(1:end-1,:);
equityEW = cumprod(ret*wEW);
equityEW = 100.*equityEW/equityEW(1);
% RP
equityRP = cumprod(ret*w_RP);
equityRP = 100.*equityRP/equityRP(1);
% Most Diversified
equityMD = cumprod(ret*w_DR);
equityMD = 100.*equityMD/equityMD(1);
% Entropy Vol Ptf
equityEVol = cumprod(ret*w_MaxEntropyVol);
equityEVol = 100.*equityEVol/equityEVol(1);

% Minimum Variance Portfolio
equityMVP = cumprod(ret*weightsMV);
equityMVP = 100.*equityMVP/equityMVP(1);

% Plot
f = figure();
plot(dates_(2:end,1), equityEW, 'LineWidth', 3)
hold on
plot(dates_(2:end,1), equityRP, 'LineWidth', 3)
hold on
plot(dates_(2:end,1), equityMD, 'LineWidth', 3)
hold on
plot(dates_(2:end,1), equityEVol, 'LineWidth', 3)
hold on
plot(dates_(2:end,1), equityMVP, 'LineWidth', 3)
legend('Equally Weighted Portfolio', 'Risk Parity Portfolio', 'Most Diversified Portfolio', 'Max Entropy Portfolio', 'Minimum Variance Portfolio')
xlabel('Date')
ylabel('Equity')

%% Performance Metrics
[annRet_ew, annVol_ew, Sharpe_ew, MaxDD_ew, Calmar_ew] = getPerformanceMetrics(equityEW);
[annRet_rp, annVol_rp, Sharpe_rp, MaxDD_rp, Calmar_rp] = getPerformanceMetrics(equityRP);
[annRet_md, annVol_md, Sharpe_md, MaxDD_md, Calmar_md] = getPerformanceMetrics(equityMD);
[annRet_ev, annVol_ev, Sharpe_ev, MaxDD_ev, Calmar_ev] = getPerformanceMetrics(equityEVol);
[annRet_MVP, annVol_MVP, Sharpe_MVP, MaxDD_MVP, Calmar_MVP] = getPerformanceMetrics(equityMVP);

metricsEW = [annRet_ew; annVol_ew; Sharpe_ew; MaxDD_ew; Calmar_ew];
metricsRP = [annRet_rp; annVol_rp; Sharpe_rp; MaxDD_rp; Calmar_rp];
metricsMD = [annRet_md; annVol_md; Sharpe_md; MaxDD_md; Calmar_md];
metricsEV = [annRet_ev; annVol_ev; Sharpe_ev; MaxDD_ev; Calmar_ev];
metricsMVP = [annRet_MVP; annVol_MVP; Sharpe_MVP; MaxDD_MVP; Calmar_MVP];

metrics_ = ["AnnRet";"AnnVola";"SharpeRatio"; "Max DD";"CalmarRatio"];
table(metrics_, metricsEW,metricsRP, metricsMD,metricsEV, metricsMVP,'VariableNames', ["Metrics","EW", "RP", "MD", "Entropy Vol", "Minimum Variance"])

%% Compare Portfolios (EW, RP, MD, EntropyVol) - OutOf Sample
start_dt = datetime('02/06/2019', 'InputFormat', 'dd/MM/yyyy'); % dt_4(1)+20;
end_dt   = datetime('01/01/2020', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed'); %fai vedere openLeft ecc 'openLeft'
subsample = myPrice_dt(rng,:);

prices_val_OS = subsample.Variables;
dates_OS = subsample.Time;
% Equity Curve
retOS = prices_val_OS(2:end,:)./prices_val_OS(1:end-1,:);
equityEW_OS = cumprod(retOS*wEW);
equityEW_OS = 100.*equityEW_OS/equityEW_OS(1);
% RP
equityRP_OS = cumprod(retOS*w_RP);
equityRP_OS = 100.*equityRP_OS/equityRP_OS(1);
% Most Diversified
equityMD_OS = cumprod(retOS*w_DR);
equityMD_OS = 100.*equityMD_OS/equityMD_OS(1);
% Entropy Vol Ptf
equityEVol_OS = cumprod(retOS*w_MaxEntropyVol);
equityEVol_OS = 100.*equityEVol_OS/equityEVol_OS(1);

% Minimum Variance Portfolio
equityMVP_OS = cumprod(retOS*weightsMV);
equityMVP_OS = 100.*equityMVP_OS/equityMVP_OS(1);

% Plot
f = figure();
plot(dates_OS(2:end,1), equityEW_OS, 'LineWidth', 3)
hold on
plot(dates_OS(2:end,1), equityRP_OS, 'LineWidth', 3)
hold on
plot(dates_OS(2:end,1), equityMD_OS, 'LineWidth', 3)
hold on
plot(dates_OS(2:end,1), equityEVol_OS, 'LineWidth', 3)
hold on
plot(dates_OS(2:end,1), equityMVP_OS, 'LineWidth', 3)
legend('Equally Weighted Portfolio', 'Risk Parity Portfolio', 'Most Diversified Portfolio', 'Max Entropy Portfolio', 'Min Variance Portfolio')
xlabel('Date')
ylabel('Equity')

%% Performance Metrics
[annRet_ew_OS, annVol_ew_OS, Sharpe_ew_OS, MaxDD_ew_OS, Calmar_ew_OS] = getPerformanceMetrics(equityEW_OS);
[annRet_rp_OS, annVol_rp_OS, Sharpe_rp_OS, MaxDD_rp_OS, Calmar_rp_OS] = getPerformanceMetrics(equityRP_OS);
[annRet_md_OS, annVol_md_OS, Sharpe_md_OS, MaxDD_md_OS, Calmar_md_OS] = getPerformanceMetrics(equityMD_OS);
[annRet_ev_OS, annVol_ev_OS, Sharpe_ev_OS, MaxDD_ev_OS, Calmar_ev_OS] = getPerformanceMetrics(equityEVol_OS);
[annRet_MVP_OS, annVol_MVP_OS, Sharpe_MVP_OS, MaxDD_MVP_OS, Calmar_MVP_OS] = getPerformanceMetrics(equityMVP_OS);

metricsEW_OS = [annRet_ew_OS; annVol_ew_OS; Sharpe_ew_OS; MaxDD_ew_OS; Calmar_ew_OS];
metricsRP_OS = [annRet_rp_OS; annVol_rp_OS; Sharpe_rp_OS; MaxDD_rp_OS; Calmar_rp_OS];
metricsMD_OS = [annRet_md_OS; annVol_md_OS; Sharpe_md_OS; MaxDD_md_OS; Calmar_md_OS];
metricsEV_OS = [annRet_ev_OS; annVol_ev_OS; Sharpe_ev_OS; MaxDD_ev_OS; Calmar_ev_OS];
metricsMVP_OS = [annRet_MVP_OS; annVol_MVP_OS; Sharpe_MVP_OS; MaxDD_MVP_OS; Calmar_MVP_OS];

metrics_OS = ["AnnRet";"AnnVola";"SharpeRatio"; "Max DD";"CalmarRatio"];
table(metrics_OS, metricsEW_OS,metricsRP_OS, metricsMD_OS,metricsEV_OS,metricsMVP_OS, 'VariableNames', ["Metrics","EW", "RP", "MD", "Entropy Vol", "Min Variance"])

clear all
close all
clc

%% Read Prices
path_map        = 'C:\Users\ginevra.angelini\Desktop\lezioni_poli\lezioni\Lezione3\';
filename        = 'geo_index_prices.xlsx';

table_prices = readtable(strcat(path_map, filename));
%% Transform prices from table to timetable
dt = table_prices(:,1).Variables;
values = table_prices(:,2:end).Variables;
nm = table_prices.Properties.VariableNames(2:end);

myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', nm); 
%% Selection of a subset of Dates
start_dt = datetime('01/05/2021', 'InputFormat', 'dd/MM/yyyy'); 
end_dt   = datetime('01/08/2021', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed'); 
subsample = myPrice_dt(rng,:);

prices_val = subsample.Variables;
dates_ = subsample.Time;
%% Calculate returns
LogRet = tick2ret(prices_val, 'Method', 'Continuous');
ExpRet = mean(LogRet);
%% Calculate Variance-Covariance Matrix
V = cov(LogRet);

%% 1.Test on Concentration Error 
p = Portfolio('AssetList', nm);
p = setDefaultConstraints(p);
P = estimateAssetMoments(p, LogRet, 'missingdata', false);
pwgt = estimateFrontier(P, 100);
[pf_Risk, pf_Retn] = estimatePortMoments(P,pwgt);
% Plot of weights
bar(pwgt', 'Stacked')
%% 1.Test on Concentration Error: Compute efficient frontier with boundaries 
LowerBound = 0.05*ones(1,8);
UpperBound = 0.8*ones(1,8);
p = setBounds(p, LowerBound, UpperBound);
MinNumAssets = 6;
MaxNumberAssets = 8;
p = setMinMaxNumAssets(p, MinNumAssets, MaxNumberAssets);
PortBound = estimateAssetMoments(p, LogRet, 'missingdata', false);
pwgt_bound = estimateFrontier(PortBound, 100);
[pf_RiskB, pf_RetnB] = estimatePortMoments(PortBound,pwgt_bound);
% Plot Weights
h = figure;
bar(pwgt', 'Stacked')
g = figure;
bar(pwgt_bound', 'Stacked')

% Plot Frontier
plot(pf_RiskB, pf_RetnB)
hold on
plot(pf_Risk, pf_Retn)
%% 2. Test the robustness of the frontier: We add small perturbations on the return of first asset
LogRet1 = LogRet;
LogRet1(:,1) = LogRet1(:,1)+std(LogRet1(:,1));

LogRet2 = LogRet;
LogRet2(:,1) = LogRet2(:,1)+1.5*std(LogRet2(:,1));

LogRet3 = LogRet;
LogRet3(:,1) = LogRet3(:,1)+2*std(LogRet3(:,1));
% Plot returns distribution
histogram(LogRet)
hold on
histogram(LogRet1, 'EdgeAlpha', 0.8)
hold on
histogram(LogRet2, 'EdgeAlpha', 0.6)
hold on 
histogram(LogRet3, 'EdgeAlpha', 0.4)

%% Test the robustness of the frontier: COMPUTE EFFICIENT FRONTIERS
p = Portfolio('AssetList', nm);
p = setDefaultConstraints(p);
Port = estimateAssetMoments(p, LogRet, 'missingdata', false);
pwgt = estimateFrontier(Port, 100);
[pf_risk, pf_Retn] = estimatePortMoments(Port, pwgt);

Port1 = estimateAssetMoments(p, LogRet1, 'missingdata', false);
pwgt1 = estimateFrontier(Port1, 100);
[pf_risk1, pf_Retn1] = estimatePortMoments(Port1, pwgt1);

Port2 = estimateAssetMoments(p, LogRet2, 'missingdata', false);
pwgt2 = estimateFrontier(Port2, 100);
[pf_risk2, pf_Retn2] = estimatePortMoments(Port2, pwgt2);

Port3 = estimateAssetMoments(p, LogRet3, 'missingdata', false);
pwgt3 = estimateFrontier(Port3, 100);
[pf_risk3, pf_Retn3] = estimatePortMoments(Port3, pwgt1);
% Plot
h = figure;
title('Portfolio Frontiers')
plot(pf_risk, pf_Retn)
hold on
plot(pf_risk1, pf_Retn1)
hold on
plot(pf_risk2, pf_Retn2)
hold on
plot(pf_risk3, pf_Retn3)
legend('Port Frontier original', 'Port Frontier 1std', 'Port Frontier 1.5std' , 'Port Frontier 2std')
xlabel('volatility')
ylabel('Expected return')
%% 3. Robust Frontier : Resampling - N simulations of returns assuming they are distibuted as a normal distribution with mean = ExpRet and covariance V
N = 100;
nAssets = 8;
RiskPtfSim = zeros(100, N);
RetPtfSim = zeros(100, N);
Weights= zeros(N, nAssets, 100); %n. sim, n.assets, n portf in a single frontier

for n = 1:N
    R = mvnrnd(ExpRet, V);
    NewExpRet = R;
    NewCov = iwishrnd(V, nAssets);
    Psim = setAssetMoments(p, NewExpRet, NewCov);
    w_sim = estimateFrontier(Psim, 100);
    [pf_riskSim, pf_RetnSim] = estimatePortMoments(Psim, w_sim);
    RetPtfSim(:,n) = pf_RetnSim;
    RiskPtfSim(:, n) = pf_riskSim;
    Weights(n,:,:) = w_sim;
end

plot(RiskPtfSim, RetPtfSim)
hold on
plot(mean(RiskPtfSim, 2), mean(RetPtfSim, 2), 'LineWidth', 4)

%% 4. Robust Frontier: Robust Estimators 
RobustExpRet = trimmean(LogRet, 10);
RobustV = robustcov(LogRet);

% Plot mean VS Trimmed Mean
x = [1,2,3,4,5,6,7,8];
m = [ExpRet; RobustExpRet];
h = figure;
bar(x,m)
legend('Original Exp Ret', 'Robust Exp Ret')

% Create Portfolio and compute frontier
p = Portfolio('AssetList', nm);
p = setDefaultConstraints(p);
Port = estimateAssetMoments(p, LogRet, 'missingdata', false);
pwgt = estimateFrontier(Port, 100);
[pf_risk, pf_Retn] = estimatePortMoments(Port, pwgt);

% Robust Mean
PortRobustMean = setAssetMoments(p, RobustExpRet, V);
pwgtRobustMean = estimateFrontier(PortRobustMean, 100);
[pf_riskRM, pf_RetnRM] = estimatePortMoments(PortRobustMean, pwgtRobustMean);

% Robust Cov
PortRobustCov = setAssetMoments(p, ExpRet, RobustV);
pwgtRobustCov = estimateFrontier(PortRobustCov, 100);
[pf_riskRC, pf_RetnRC] = estimatePortMoments(PortRobustCov, pwgtRobustCov);

% Robust Cov & Robust Mean
PortRobust = setAssetMoments(p, RobustExpRet, RobustV);
pwgtRobust = estimateFrontier(PortRobust, 100);
[pf_riskRobust, pf_RetnRobust] = estimatePortMoments(PortRobust, pwgtRobust);

% Plot
h = figure;
plot(pf_risk, pf_Retn)
hold on
plot(pf_riskRM, pf_RetnRM)
hold on
plot(pf_riskRC, pf_RetnRC)
hold on
plot(pf_riskRobust, pf_RetnRobust)
legend('Port Frontier', 'PF-Robust Mean', 'PF- Robust Cov', 'Robust Port Frontier')
%% 5. Robust Frontier : Shrinkage Estimators -> Bayes-Stein Estimator
N = 8;
T = length(LogRet);
e = ExpRet;
lambda = ((N+2)*(T-1))/((e-ones(1,N).*mean(e))*inv(V)*(e-ones(1,N).*mean(e))'.*(T-N-2));
alpha = lambda/(lambda+T);
BSmean = (1-alpha).*e+alpha.*ones(1,N).*mean(e);
% all weights sum to 1, no shorting, and 100% investment in risky assets
p = Portfolio('AssetList', nm);
p = setDefaultConstraints(p);
PortBS = setAssetMoments(p, BSmean, V);
pwgtBS = estimateFrontier(PortBS, 100);
[pf_riskBS, pf_RetnBS] = estimatePortMoments(PortBS,pwgtBS );

%% 6. Robust Frontier : Compare Frontiers
h = figure;
plot(pf_risk, pf_Retn, 'LineWidth', 2)
hold on
plot(pf_riskBS, pf_RetnBS, 'LineWidth', 2)
hold on
plot(pf_riskRobust, pf_RetnRobust, 'LineWidth', 2)
hold on
plot(mean(pf_riskSim,2), mean(pf_RetnSim,2), 'LineWidth', 2)
legend('Original Frontier', 'bayes-stein Frotier', 'Robust Estimator Frontier', 'resampled Frontier')

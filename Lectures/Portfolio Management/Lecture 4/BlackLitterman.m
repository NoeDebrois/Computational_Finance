clear all
close all
clc

%% Read Prices
load myPrice_dt
load array_prices

%% Calculate returns and Covariance Matrix
Ret = tick2ret(prices_val);
numAssets = size(Ret, 2);
CovMatrix = cov(Ret);
%% Building the views
v = 3; % total 3 views v = num views
tau = 1/length(Ret);
P = zeros(v, numAssets);
q = zeros(v, 1);
Omega = zeros(v);
% View 1: 5% annual return of Apple
P(1, assetNames == 'APPLUSEquity') =1;
q(1) = 0.05;
% View 2: 3% annual return of Amazon
P(2, assetNames == 'AMZNUSEquity') =1;
q(2) = 0.03;
% View 3 Google is going to outperform JPMorgan by 5% 
P(3, assetNames == 'GOOGLUSEquity') = 1;
P(3, assetNames == 'JPMUSEquity') = -1;
q(3) = 0.05;

Omega(1,1) = tau.*P(1,:)*CovMatrix*P(1,:)';
Omega(2,2) = tau.*P(2,:)*CovMatrix*P(2,:)';
Omega(3,3) = tau.*P(3,:)*CovMatrix*P(3,:)';
% from annual view to daily view
bizyear2bizday = 1/252;
q = q*bizyear2bizday;
Omega = Omega*bizyear2bizday;

% Plot views distribution
X_views = mvnrnd(q, Omega, 200);
histogram(X_views)
%% market implied ret
% the equilibrium returns are likely equal to the implied returns from the equilibrium portfolio holding. 
%In practice, the applicable equilibrium portfolio holding can be any optimal portfolio that the investment analyst would use 
%in the absence of additional views on the market, such as the portfolio benchmark, an index, or even the current portfolio
load cap
wMKT = cap(1:15)/sum(cap(1:15));
lambda = 1.2;
mu_mkt = lambda.*CovMatrix*wMKT;
C = tau.*CovMatrix;
% plot prior distribution
X = mvnrnd(mu_mkt, C, 200);
histogram(X)
%% Black Litterman
muBL = inv(inv(C)+P'*inv(Omega)*P)*(P'*inv(Omega)*q + inv(C)*mu_mkt); 
covBL = inv(P'*inv(Omega)*P + inv(C));

table(assetNames', mu_mkt*252, muBL*252, 'VariableNames', ["AssetNames", "Prior Belief on Exp Ret", "BL ExpREt"])
% Plot Distribution

%% Portfolio Optimization
port = Portfolio('NumAssets', numAssets, 'Name', 'MeanVariance');
port = setDefaultConstraints(port);
Port = setAssetMoments(port, mean(Ret), CovMatrix);
pwgt = estimateFrontier(Port, 100);
[pf_vola, pf_ret] = estimatePortMoments(Port, pwgt);
rf = 0;
% MaxSharpe
N = 100;
ExpRet = mean(Ret);
SharpeArray = ones(1,N);
for i = 1:N
    s = getSharpeRatio(pwgt(:,i), ExpRet, CovMatrix, rf);
    SharpeArray(i) = s;
end

[MaxSharpe, Index_] = max(SharpeArray);
w_maxsharpe = pwgt(:, Index_);
wts = estimateMaxSharpeRatio(Port);

sum(wts-w_maxsharpe)
% Plot Sharpe vs vola
plot(pf_vola, SharpeArray, '-o', 'LineWidth', 3)
hold on
scatter (pf_vola(Index_), MaxSharpe, 'LineWidth', 3)


%% Black-Litterman PTF
portBL = Portfolio('NumAssets', numAssets, 'Name', 'MV with BL');
portBL = setDefaultConstraints(portBL);
portBL = setAssetMoments(portBL, muBL, CovMatrix+covBL);
wtsBL = estimateMaxSharpeRatio(portBL);
[risk, ret] = estimatePortMoments(portBL, wtsBL);

% Plot
ax1 = subplot(1,2,1);
idx = wts >0.001;
pie(ax1, wts(idx), assetNames(idx));
title(ax1, port.Name, 'Position', [-0.05, 1.6, 0]);

ax2 = subplot(1,2,2);
idx_BL = wtsBL >0.001;
pie(ax2, wtsBL(idx_BL), assetNames(idx_BL));
title(ax2, portBL.Name, 'Position', [-0.05, 1.6, 0]);

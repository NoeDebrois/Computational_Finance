clear all
close all
clc

%% Read Prices
path_map        = 'C:\Users\ginevra.angelini\Desktop\lezioni_poli\lezioni\Lezione2\';
filename        = 'asset_prices_student.xlsx';

table_prices = readtable(strcat(path_map, filename));
%% Transform prices from table to timetable
dt = table_prices(:,1).Variables;
values = table_prices(:,2:end).Variables;
nm = table_prices.Properties.VariableNames(2:end);

myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', nm); %già in datetime
%% Selection of a subset of Dates
start_dt = datetime('01/02/2021', 'InputFormat', 'dd/MM/yyyy'); % dt_4(1)+20;
end_dt   = datetime('01/04/2021', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed'); %fai vedere openLeft ecc 'openLeft'
subsample = myPrice_dt(rng,:);

prices_val = subsample.Variables;
dates_ = subsample.Time;
%save myPrice_dt myPrice_dt
%% Calculate returns
% Method 1 
ret = prices_val(2:end,:)./prices_val(1:end-1,:);
LogRet = log(ret);

% method 2 
%LogRet1 = tick2ret(prices_val, 'Method','Continuous');

ExpRet = mean(LogRet);
%% Calculate Variance
var_ = var(LogRet);
std_ = std(LogRet);

V = cov(LogRet);
%% Creation of N random portfolio
N = 100000;
NumAssets = 15;
RetPtfs = zeros(1,N);
VolaPtfs = zeros(1,N);
SharpePtfs = zeros(1,N);

for n = 1:N
    w = rand(1,NumAssets);
    w_norm = w./sum(w);
    
    exp_ret_ptf = w_norm*ExpRet';
    exp_vola_ptf = sqrt(w_norm*V*w_norm');
    sharpe_ratio = exp_ret_ptf/exp_vola_ptf;

    RetPtfs(n) = exp_ret_ptf;
    VolaPtfs(n) = exp_vola_ptf;
    SharpePtfs(n) = sharpe_ratio;
end

%% Plot
h=figure;
title('Expected return vs volatility')
scatter(VolaPtfs, RetPtfs, [], SharpePtfs, 'filled')
colorbar
xlabel('Volatility')
ylabel('Expected return')


%% Portfolio Frontier
fun = @(x)x'*V*x;
ret_ = linspace(min(RetPtfs), max(RetPtfs),100);
%ret_ = linspace(0,0.005,100)
x0 = rand(1,NumAssets)';
x0 = x0/sum(x0);
lb = zeros(1,NumAssets);
ub = ones(1,NumAssets);
FrontierVola = zeros(1, length(ret_));
FrontierRet = zeros(1, length(ret_));

for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r];
    %find best vola
    w_opt = fmincon(fun, x0, [], [], Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt'*V*w_opt);
    
    FrontierVola(i) = min_vola;
    FrontierRet(i)= r; %w_opt'*exp_ret';
end
%% Plot
h=figure;
title('Expected return vs volatility')
scatter(VolaPtfs, RetPtfs, [], SharpePtfs, 'filled')
colorbar
hold on
plot(FrontierVola, FrontierRet)
xlabel('Volatility')
ylabel('Expected return')

%% Add Costraints: min exp = 0.01 and max exp = 0.7 for each asset 
fun = @(x)x'*V*x;
ret_ = linspace(min(RetPtfs), max(RetPtfs),100);
x0 = rand(NumAssets, 1);
x0 = x0/sum(x0);
lb = zeros(1,NumAssets);
ub = ones(1,NumAssets);

% Second Costraint
A_max = eye(NumAssets);
A_min= -eye(NumAssets);
b_min = -0.01.*ones(NumAssets,1);
b_max = 0.7*ones(NumAssets,1);
A = [A_min; A_max];
b = [b_min, b_max];

FrontierVola2 = zeros(1, length(ret_));
FrontierRet2 = zeros(1, length(ret_));

for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r];
    % find optimal w, minimizing volatility
    w_opt= fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt'*V*w_opt);
    
    FrontierVola2(i) = min_vola;
    FrontierRet2(i)= r; %w_opt'*exp_ret';
end

%% Add Costraints - Part 2 : AAPL, AMZN and GOOGL exposition has to be more than 10%
fun = @(x)x'*V*x;
ret_ = linspace(min(RetPtfs), max(RetPtfs),100);
x0 = rand(NumAssets, 1);
x0 = x0./sum(x0);
lb = zeros(1,NumAssets);
ub = ones(1,NumAssets);

% Second Costraint
A = [-1,0,-1,0,0,0,-1,0,0,0,0,0,0,0,0];
b = -0.25;

FrontierVola3 = zeros(1, length(ret_));
FrontierRet3 = zeros(1, length(ret_));
d = zeros(1,length(ret_));
for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r];
    % find optimal w, minimizing volatility
    w_opt= fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt'*V*w_opt);
    
    FrontierVola3(i) = min_vola;
    FrontierRet3(i)= r; %w_opt'*exp_ret';
    d(i) =w_opt(1)+w_opt(3)+w_opt(7);
end

%% Plot the three frontiers
WeightsEW = 1/15.*ones(1,15);
h=figure;
title('Expected return vs volatility')
plot(FrontierVola, FrontierRet,'LineWidth', 4)
hold on
plot(FrontierVola2, FrontierRet2, 'LineWidth',3)
hold on
plot(FrontierVola3, FrontierRet3,'LineWidth',3)
hold on
scatter(sqrt(WeightsEW*V*WeightsEW'), WeightsEW*ExpRet', 'filled')
legend('Frontier 1', 'Frontier 2', 'Frontier 3', 'EW Ptf')
xlabel('Volatility')
ylabel('Expected return')

%% Portfolio Frontier respect to benchmark
WeightsEW = 1/NumAssets.*ones(1,NumAssets);
VolaEW = sqrt(WeightsEW*V*WeightsEW');
RetEW = WeightsEW*ExpRet';
fun = @(x)(x'-WeightsEW)*V*(x-WeightsEW');
ret_ = linspace(RetEW, max(RetPtfs)*2,100);
x0 = rand(1,NumAssets)';
x0 = x0/sum(x0);
lb = zeros(1,NumAssets);
ub = ones(1,NumAssets);
FrontierVolaBench = zeros(1, length(ret_));
FrontierRetBench = zeros(1, length(ret_));
Weights = zeros(NumAssets, length(ret_));
NumAssetInPort = zeros(1, length(ret_));

for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r+(WeightsEW*ExpRet')];
    %find best vola
    [w_opt, fval] = fmincon(fun, x0, [], [], Aeq, beq, lb, ub);
    min_vola = sqrt(fval); 
    
    FrontierVolaBench(i) = min_vola;
    FrontierRetBench(i)= r;
    Weights(:,i) = w_opt;
    NumAssetInPort(i) = length(w_opt(round(w_opt,2) ~= 0));
end

% Calculus of relative return and relative volatility
FrontierRetRel = (FrontierRetBench./RetEW);
FrontierVolaRel = (FrontierVolaBench./VolaEW);

% Plot
scatter(FrontierVolaRel, FrontierRetRel, 'filled', 'LineWidth', 1)
ylim([0, max(FrontierRetRel)+0.5])
xlabel('Relative Volatility')
ylabel('Relative Expected Return')

% Plot 
f = figure();
histogram(NumAssetInPort)
%% Portfolio object
p = Portfolio('AssetList',nm);

% all weights sum to 1, no shorting, and 100% investment in risky assets).
p = setDefaultConstraints(p);
P = estimateAssetMoments(p, LogRet,'missingdata',false);

%% COMPUTE EFFICIENT FRONTIER AND PLOT INFORMATION RATIO (IRR)
pwgt = estimateFrontier(P, 100); % Estimate the weights.
[pf_Risk, pf_Retn] = estimatePortMoments(P, pwgt); % Get the risk and return.

%% Add turnover costraints and Transaction costs
BuyCost = 0.0020;
SellCost = 0.0020;
Turnover = 0.3; %average turnover < 30%

p = setInitPort(p,WeightsEW);
q = setCosts(p, BuyCost,SellCost);
q = setTurnover(q,Turnover);
P1 = estimateAssetMoments(q,LogRet,'missingdata',false);
pwgt1 = estimateFrontier(P1, 100);
[pf_Risk1, pf_Retn1] = estimatePortMoments(P1, pwgt1); % Get the risk and return.

% Plot
h = figure;
title('Expected Return VS Volatility')
plot(pf_Risk,pf_Retn, 'LineStyle', '--', 'LineWidth', 2)
hold on
plot(pf_Risk1, pf_Retn1, 'LineWidth', 2)
xlabel('Volatility')
ylabel('Expected Return')
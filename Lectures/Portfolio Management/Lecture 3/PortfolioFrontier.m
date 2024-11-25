% Portfolio Management - Lesson 3.
%% Portfolio frontiers with or without constraints
clear
close
clc
%
%% Read Prices :
path_map        = '/Users/noedebrois/Desktop/Desktop - Noé’s MacBook Air/Politecnico/Computational Finance/AY 2024:2025/Portfolio Management/Lesson 3/';
filename        = 'asset_prices_lesson2.xlsx';
table_prices = readtable(strcat(path_map, filename));
%
%% Transform prices from table to timetable :
dt = table_prices(:,1).Variables; % All the dates
values = table_prices(:,2:end).Variables; % All the prices (starting at column 2)
nm = table_prices.Properties.VariableNames(2:end); % Column names
%
% Translate to timetable :
% - Reminder : timetables = tables for time series data, with timestamped 
% rows and variables of different types.
myPrice_dt = array2timetable(values, 'RowTimes', dt, 'VariableNames', nm);
%
%% Selection of a subset of Dates :
start_dt = datetime('01/02/2021', 'InputFormat', 'dd/MM/yyyy');
end_dt   = datetime('01/04/2021', 'InputFormat', 'dd/MM/yyyy');
% 
rng = timerange(start_dt, end_dt, 'closed'); % Look at rng to see what we selected
subsample = myPrice_dt(rng,:); % Still a timetable !
% 
prices_val = subsample.Variables; % Get the prices from the subsample
dates_ = subsample.Time; % Get the time from the subsample
%
%% Calculate daily log-returns :
% Method 1 :
ret = prices_val(2:end,:)./prices_val(1:end-1,:); % Returns
LogRet = log(ret); % Log-returns
%
% Method 2 : with 'tick2ret' function :
% - Help : "Convert price series to return series"
LogRet1 = tick2ret(prices_val, 'Method','Continuous'); % It gives the same!
% Compute the average daily log-return for each equity stock :
ExpRet = mean(LogRet);
%
%% Calculate Variance :
var_ = var(LogRet); % Variance
std_ = std(LogRet); % Standard-deviation
%
V = cov(LogRet); % Covariance matrix 
%
%% Creation of N random portfolios :
N = 100000;
NumAssets = 15;
% Arrays to store the results of our generated portfolios :
RetPtfs = zeros(1,N); % Expected Returns
VolaPtfs = zeros(1,N); % Volatilities
SharpePtfs = zeros(1,N); % Sharpe ratios
%
for n = 1:N
    w = rand(1,NumAssets); % Random ptf-weights
    w_norm = w./sum(w); % Normalize ptf-weights
    % 
    exp_ret_ptf = w_norm * ExpRet'; % Expected return of this random ptf
    exp_vola_ptf = sqrt(w_norm * V * w_norm'); % Volatility of this random ptf
    sharpe_ratio = exp_ret_ptf / exp_vola_ptf; % Sharpe ratio of this random ptf
    % Put everything into arrays :
    RetPtfs(n) = exp_ret_ptf;
    VolaPtfs(n) = exp_vola_ptf;
    SharpePtfs(n) = sharpe_ratio;
end
%
%% Plot :
% WE PLOT JUST AFTER THE NEXT BLOCK
% figure;
% scatter(VolaPtfs, RetPtfs, [], SharpePtfs, 'filled')
% colorbar
% xlabel('Volatility')
% ylabel('Expected return')
% title('Expected return vs volatility [for our N random portfolios]')
%
%% Portfolio Frontier :
fun = @(x)x' * V * x; % Function defined directly here
ret_ = linspace(min(RetPtfs),max(RetPtfs),100);
%
x0 = rand(1,NumAssets)'; % Uniformy-distributed random vector 1 by NumAssets
x0 = x0 / sum(x0); % Normalise
lb = zeros(1,NumAssets); % Lower bound
ub = ones(1,NumAssets); % Upper bound
FrontierVola = zeros(1, length(ret_));
FrontierRet = zeros(1, length(ret_));
%
for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq = [1; r];
    % Find best volatility :
    % - Help : fmincon "finds minimum of constrained nonlinear 
    % multivariable function"
    w_opt = fmincon(fun, x0, [], [], Aeq, beq, lb, ub); % Minimisation algo
    min_vola = sqrt(w_opt' * V * w_opt); % Compute the vol w/ the optimum weights
    % 
    FrontierVola(i) = min_vola; % For this return, add the minimum volatility
    FrontierRet(i) = r; % Add the return ; w_opt' * exp_ret';
end
%
%% Plot :
figure;
scatter(VolaPtfs, RetPtfs, [], SharpePtfs, 'filled')
colorbar
hold on
plot(FrontierVola, FrontierRet)
xlabel('Volatility')
ylabel('Expected return')
title('Expected return vs volatility [for our N random portfolios] & the portfolio frontier')
%
%% Add Constraints - PART 1 :
% min exp = 0.01 and max exp = 0.7 for each asset 
fun = @(x)x' * V * x; % Same as above
ret_ = linspace(min(RetPtfs), max(RetPtfs),100); % Same as above
x0 = rand(NumAssets, 1); % Same as above
x0 = x0/sum(x0); % Same as above
lb = zeros(1,NumAssets); % Same as above
ub = ones(1,NumAssets); % Same as above
% Constraints :
A_max = eye(NumAssets);
A_min = -eye(NumAssets);
b_min = -0.01.*ones(NumAssets,1);
b_max = 0.7*ones(NumAssets,1);
A = [A_min; A_max];
b = [b_min, b_max];
%
FrontierVola2 = zeros(1, length(ret_));
FrontierRet2 = zeros(1, length(ret_));
% Same as previously, but with the two constraints (inside A and b):
for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r];
    % Find optimal w, minimizing volatility
    w_opt = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt'*V*w_opt);
    % Update the frontier :
    FrontierVola2(i) = min_vola;
    FrontierRet2(i)= r; %w_opt'*exp_ret';
end
%
%% Add Constraints - Part 2 : 
% AAPL, AMZN and GOOGL exposition has to be more than 10%
% Same as before :
fun = @(x)x' * V * x;
ret_ = linspace(min(RetPtfs), max(RetPtfs),100);
x0 = rand(NumAssets,1);
x0 = x0./sum(x0);
lb = zeros(1,NumAssets);
ub = ones(1,NumAssets);
% Second constraint :
A = [-1,0,-1,0,0,0,-1,0,0,0,0,0,0,0,0];
b = -0.25;
% Same as before but with the added constraint (inside A and b) :
FrontierVola3 = zeros(1, length(ret_));
FrontierRet3 = zeros(1, length(ret_));
d = zeros(1,length(ret_));
for i = 1:length(ret_)
    r = ret_(i);
    Aeq = [ones(1,NumAssets); ExpRet]; 
    beq =[1; r];
    % Find optimal w, minimizing volatility
    w_opt = fmincon(fun, x0, A, b, Aeq, beq, lb, ub);
    min_vola = sqrt(w_opt' * V * w_opt);
    % Update the frontier :
    FrontierVola3(i) = min_vola;
    FrontierRet3(i)= r; %w_opt'*exp_ret';
    d(i) = w_opt(1) + w_opt(3) + w_opt(7);
end
%
%% Plot the three frontiers :
WeightsEW = 1/15.*ones(1,15); % Equilibrated portfolio (flat, i.e : all the same weights for each asset)
figure;
plot(FrontierVola, FrontierRet,'LineWidth', 4)
hold on
plot(FrontierVola2, FrontierRet2, 'LineWidth',3)
hold on
plot(FrontierVola3, FrontierRet3,'LineWidth',3)
hold on
scatter(sqrt(WeightsEW*V*WeightsEW'), WeightsEW*ExpRet', 'filled')
legend('Frontier 1 - no constraint', 'Frontier 2 - constraint on min & max weight', 'Frontier 3 - constrain on AAPL, AMZN and GOOGL weights', 'EW Ptf')
xlabel('Volatility')
ylabel('Expected return')
title('Expected return vs volatility - Portfolio frontiers with different constraints')
%
%% Portfolio Frontier with respect to benchmark (flat weights benchmark)
%
% Step 1: Set up the benchmark portfolio (equal weights across assets)
WeightsEW = 1/NumAssets.*ones(1,NumAssets);  % Equal weights for each asset
VolaEW = sqrt(WeightsEW * V * WeightsEW');   % Volatility of the benchmark portfolio
RetEW = WeightsEW * ExpRet';                 % Expected return of the benchmark portfolio
%
% Step 2: Define the objective function for volatility minimization
fun = @(x)(x' - WeightsEW) * V * (x - WeightsEW');  % Function to minimize volatility relative to the benchmark
%
% Step 3: Set up the range of returns for which we want to find the optimal portfolio
ret_ = linspace(RetEW, max(RetPtfs)*2, 100);  % Returns range from the benchmark to twice the maximum return
x0 = rand(1, NumAssets)';                     % Initial random guess for asset weights
x0 = x0 / sum(x0);                            % Normalize the weights so that their sum equals 1
lb = zeros(1, NumAssets);                     % Lower bound (no short selling, weights ≥ 0)
ub = ones(1, NumAssets);                      % Upper bound (weights ≤ 1)
%
% Step 4: Initialize arrays to store results
FrontierVolaBench = zeros(1, length(ret_));   % Array to store the volatilities of the optimal portfolios
FrontierRetBench = zeros(1, length(ret_));    % Array to store the expected returns of the optimal portfolios
Weights = zeros(NumAssets, length(ret_));     % Array to store the optimal weights of the assets
NumAssetInPort = zeros(1, length(ret_));      % Array to store the number of significant assets in each portfolio
%
% Step 5: Loop over each target return, find the optimal portfolio
for i = 1:length(ret_)
    r = ret_(i);  % Target return for this iteration
    Aeq = [ones(1, NumAssets); ExpRet];   % Equality constraints: sum(weights) = 1 and expected return = r
    beq = [1; r + (WeightsEW * ExpRet')]; % Adjust the target return with respect to the benchmark
    %
    % Find the portfolio with minimum volatility for the target return
    [w_opt, fval] = fmincon(fun, x0, [], [], Aeq, beq, lb, ub); 
    min_vola = sqrt(fval);  % The minimum volatility (sqrt of the function value)
    %
    % Store the results for this iteration
    FrontierVolaBench(i) = min_vola;  % Store the optimal volatility
    FrontierRetBench(i) = r;          % Store the corresponding return
    Weights(:,i) = w_opt;             % Store the optimal asset weights
    NumAssetInPort(i) = length(w_opt(round(w_opt,2) ~= 0));  % Count the number of non-zero weights
end
%
% Step 6: Calculate relative return and relative volatility compared to the benchmark
FrontierRetRel = (FrontierRetBench ./ RetEW);     % Relative return (compared to benchmark return)
FrontierVolaRel = (FrontierVolaBench ./ VolaEW);  % Relative volatility (compared to benchmark volatility)
%
% Step 7: Plot relative return vs. relative volatility
figure;
scatter(FrontierVolaRel, FrontierRetRel, 'filled', 'LineWidth', 1)  % Scatter plot
ylim([0, max(FrontierRetRel) + 0.5])  % Set y-axis limits
xlabel('Relative Volatility')         % Label for x-axis
ylabel('Relative Expected Return')    % Label for y-axis
title("Relative return vs relative volatility for the benchmark allocation (flat weights)")
hold off
%
% Step 8: Plot the histogram of the number of assets in the optimal portfolios
f = figure();                        % Create a new figure
histogram(NumAssetInPort)            % Plot histogram of the number of assets in each portfolio
title("Histogram of the number of assets in each portfolio")
hold off
%
%% Portfolio object creation:
p = Portfolio('AssetList', nm);  % Create a Portfolio object with a list of assets ('nm' is the list of asset names)

% Set default constraints: sum of weights equals 1 (fully invested), no shorting allowed
p = setDefaultConstraints(p);

% Estimate asset moments (mean returns and covariance) from the log returns data 'LogRet'
% 'missingdata', false ensures no handling of missing data (assumes no missing values)
P = estimateAssetMoments(p, LogRet, 'missingdata', false);
%
%% COMPUTE EFFICIENT FRONTIER AND PLOT INFORMATION RATIO (IRR):

% Estimate the efficient frontier weights for 100 points along the frontier
pwgt = estimateFrontier(P, 100);

% Get the portfolio risk (volatility) and expected return for each point on the frontier
[pf_Risk, pf_Retn] = estimatePortMoments(P, pwgt);
%
%% Add turnover constraints and transaction costs:

% Define transaction costs for buying and selling assets
BuyCost = 0.0020;  % 0.20% cost for buying
SellCost = 0.0020; % 0.20% cost for selling

% Set a turnover constraint: average turnover must be less than 30%
Turnover = 0.3;  % 30% turnover limit

% Set the initial portfolio to the benchmark portfolio with equal weights (WeightsEW)
p = setInitPort(p, WeightsEW);

% Add the transaction costs for buying and selling
q = setCosts(p, BuyCost, SellCost);

% Set the turnover constraint in the portfolio object
q = setTurnover(q, Turnover);

% Estimate asset moments (mean returns and covariance) again, with transaction costs and turnover constraint
P1 = estimateAssetMoments(q, LogRet, 'missingdata', false);

% Estimate the efficient frontier under the new constraints (100 points along the frontier)
pwgt1 = estimateFrontier(P1, 100);

% Get the portfolio risk (volatility) and expected return for the new frontier
[pf_Risk1, pf_Retn1] = estimatePortMoments(P1, pwgt1);
%
%% Plot the efficient frontier before and after adding constraints:

h = figure;  % Create a new figure for the plot

% Plot the original efficient frontier (without constraints)
plot(pf_Risk, pf_Retn, 'LineStyle', '--', 'LineWidth', 2)  % Dashed line

hold on  % Keep the plot open to overlay the new frontier

% Plot the efficient frontier with turnover constraints and transaction costs
plot(pf_Risk1, pf_Retn1, 'LineWidth', 2)  % Solid line

% Legend :
legend('Original efficient frontier (w/o constraint)', 'Efficient frontier w/ turnover constraint & transaction costs')

% Add labels to the axes and title
xlabel('Volatility')  % Label for the x-axis (Risk)
ylabel('Expected Return')  % Label for the y-axis (Return)
title('Expected Return VS Volatility')  % Set the title of the plot

% We see that, for a given return, the variance of the original efficient
% frontier (with no constraint) is smaller than the variance of the
% efficient frontier with constraints. This is normal, with constraint we
% get a "suboptimal" solution, but on the other hand we can reach some
% goals (like diversification for example).
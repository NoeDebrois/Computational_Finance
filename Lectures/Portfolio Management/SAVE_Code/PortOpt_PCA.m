clear all
close all
clc

%% Read Prices & Names
load array_prices.mat
load myPrice_dt.mat
start_dt = datetime('01/01/2021', 'InputFormat', 'dd/MM/yyyy'); % dt_4(1)+20;
end_dt   = datetime('01/01/2022', 'InputFormat', 'dd/MM/yyyy');

rng = timerange(start_dt, end_dt,'closed'); %fai vedere openLeft ecc 'openLeft'
subsample = myPrice_dt(rng,:);

prices_val = subsample.Variables;
dates_ = subsample.Time;
%% Compute Moments
LogRet = tick2ret(prices_val, 'Method','Continuous');
ExpRet = mean(LogRet);
RetStd = (LogRet-ExpRet)./ std(LogRet);

CovMatrix = cov(LogRet);
%xlswrite('C:\Users\ginevra.angelini\Desktop\poli_secondo_anno\lezioni\Lezione7\log_ret.xlsx',LogRet)

%% pca
k = 10;
[factorLoading,factorRetn,latent,r, explained, mu] = pca(RetStd, 'NumComponents', k);
covarFactor = cov(factorRetn);

TotVar = sum(latent);
ExplainedVar = latent(1:k)/TotVar;
n_list = [1,2,3,4,5,6,7,8,9,10];
CumExplainedVar = zeros(1, size(n_list,2));
for i =1:size(n_list,2)
    n = n_list(i);
    CumExplainedVar(1,i) = getCumulativeExplainedVar_old(latent, n);
end

CumExplVar_copy = cumsum(explained);
%CumExplVar_copy(1:10)'
%CumExplainedVar*100

% Plot
h = figure();
title('Percentage of Explained Variances for each Principal Component')
bar(n_list,ExplainedVar)
xlabel('Principal Components')
ylabel('Percentage of Explained Variances')

% plot 2
f = figure();
title('Total Percentage of Explained Variances for the first n-components')
plot(n_list,CumExplainedVar, 'm')
hold on
scatter(n_list,CumExplainedVar,'m', 'filled')
grid on
xlabel('Total number of Principal Components')
ylabel('Percentage of Explained Variances')

%reconstruct asset returns
reconReturn = factorRetn*factorLoading' + ExpRet;
unexplainedRetn = LogRet - reconReturn; % epsilon

% There are unexplained asset returns εa because the remaining (p - k) principal components are dropped. 
% You can attribute the unexplained asset returns to the asset-specific risks represented as D.
unexplainedCovar = diag(cov(unexplainedRetn));
D = diag(unexplainedCovar);


covarAsset = factorLoading*covarFactor*factorLoading'+D;
%% Optimization max(ret-variance)
func = @(x) -((ExpRet*x)-((factorLoading'*x)'*covarFactor*(factorLoading'*x)+x'*D*x));

x0 = rand(size(LogRet,2),1);
x0 = x0./sum(x0);
lb = zeros(1,size(LogRet,2));
ub = ones(1,size(LogRet,2));
Aeq = ones(1,size(LogRet,2));
beq = 1;
[w_opt, fval] = fmincon(func, x0, [],[],Aeq, beq, lb, ub);
wf = (factorLoading'*w_opt);

% equity curve & performances
ret = prices_val(2:end,:)./prices_val(1:end-1,:);
equity_ = cumprod(ret*w_opt);
equity_ = 100.*equity_/equity_(1);
[annRet_, annVol_, Sharpe_, MaxDD_, Calmar_] = getPerformanceMetrics(equity_);
perfTable_ = table(annRet_, annVol_, Sharpe_, MaxDD_, Calmar_, 'VariableNames',["AnnRet", "AnnVol", "Sharpe", "MaxDD","Calmar"]);

% Plot Equity pca VS equity EW
wEW = 1/15*ones(15, 1);
equity_ew = cumprod(ret*wEW);
equity_ew = 100.*equity_ew/equity_ew(1);
[annRet_ew, annVol_ew, Sharpe_ew, MaxDD_ew, Calmar_ew] = getPerformanceMetrics(equity_ew);
perfTable_ew = table(annRet_ew, annVol_ew, Sharpe_ew, MaxDD_ew, Calmar_ew, 'VariableNames',["AnnRet", "AnnVol", "Sharpe", "MaxDD","Calmar"]);

% Plot
f1 = figure();
plot(dates_(2:end,1),equity_)
hold on
plot(dates_(2:end,1),equity_ew)
legend('Equity PCA', 'Equity EW')
%% Sensitivity Analysis + Optimization MaxSharpe
k_list = [2,5,7,10,12];
x0 = rand(size(LogRet,2),1);
x0 = x0./sum(x0);
lb = zeros(1,size(LogRet,2));
ub = ones(1,size(LogRet,2));
Aeq = ones(1,size(LogRet,2));
beq = 1;
weights_pca = zeros(size(LogRet,2), size(k_list,2));
equity_matrix = zeros(size(equity_ew,1), size(k_list,2));
ExplVar = ones(size(k_list,2),1);
perf_table =  table('Size',[5 5],'VariableTypes',{'double', 'double','double','double','double'}, 'VariableNames',{'AnnRet', 'AnnVol', 'Sharpe', 'MaxDD','Calmar'},'RowNames',{'k = 2','k = 5','k = 7', 'k = 10','k = 12'});

for i= 1:size(k_list,2)
    k = k_list(i);
    [factorLoading,factorRetn,latent,~, explained] = pca(RetStd, 'NumComponents', k);
    cumulativeExplained = cumsum(explained);
    ExplVar(i,:) = cumulativeExplained(k);
    covarFactor = cov(factorRetn);
    reconReturn = factorRetn*factorLoading' + ExpRet;
    unexplainedRetn = LogRet - reconReturn; % epsilon
    unexplainedCovar = diag(cov(unexplainedRetn));
    D = diag(unexplainedCovar);
    covarAsset = factorLoading*covarFactor*factorLoading'+D;
    f_max_sharpe = @(x) -((ExpRet*x)/sqrt((factorLoading'*x)'*covarFactor*(factorLoading'*x)+x'*D*x));
    [w_opt, fval] = fmincon(f_max_sharpe, x0, [],[],Aeq, beq, lb, ub);
    weights_pca(:,i) = w_opt;
    equity_temp = cumprod(ret*w_opt);
    equity_matrix(:,i)= 100.*equity_temp/equity_temp(1);
    [annRet_, annVol_, Sharpe_, MaxDD_, Calmar_] = getPerformanceMetrics(equity_temp);
    perf_table(i,:) = table(annRet_, annVol_, Sharpe_, MaxDD_, Calmar_);
end

%plot Expl Var
f1 = figure();
plot(k_list, ExplVar)
hold on 
scatter(k_list, ExplVar, 'filled')

% Plot Equity
f2 = figure();
plot(dates_(2:end,1), equity_matrix)
hold on 
plot(dates_(2:end,1),equity_ew, 'k')
grid on
legend('k = 2','k = 5', 'k = 7', 'k = 10','k =12', 'EW')
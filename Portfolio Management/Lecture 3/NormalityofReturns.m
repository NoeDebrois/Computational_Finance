% Portfolio Management - Lesson 3.
%% Are returns normal ? This file helps to understand that the answer is NO
%
clear all
close all
clc
%
%% Read Prices :
path_map        = '/Users/noedebrois/Desktop/Desktop - Noé’s MacBook Air/Politecnico/Computational Finance/AY 2024:2025/Portfolio Management/Lesson 3/';
filename        = 'geo_index_prices.xlsx';
table_prices = readtable(strcat(path_map, filename));
%
%% Calculate Returns (all period) :
dt = table_prices(:,1).Variables; % All the dates from the file
values = table_prices(:,2).Variables; % The second column (DAX index)
% WARNING : the first 2 values are NAN, hence the start at 4 and 3 instead
% of 2 and 1 :
Ret = values(4:end)./values(3:end-1); % Compute the daily return P_end_of_the_day/P_beginning_of_the_day
% Plot
h = figure;
histogram(Ret, 'EdgeAlpha', 0.8)
title("Daily returns of the DAX index")
%
%% Generation of N  RV ~ normal distribution with same mean and variance :
N = length(Ret);
NormDistrVar = zeros(1,N);
for i = 1:N
    r_ = normrnd(mean(Ret),std(Ret));
    NormDistrVar(i) = r_;
end
%
%% Ensemble plot :
h = figure;
histogram(Ret, 'EdgeAlpha', 0.8)
hold on 
histogram(NormDistrVar, 'EdgeAlpha', 0.7)
legend('Original Returtns Distr', 'Generated Values')
title("Comparison between empirical daily returns & normally-distributed returns (with same mean & variance)")
%
% CLEARLY WE SEE THAT EXTREME VALUES HAPPEN IN HISTORICAL DATA, BUT DO NOT
% HAPPEN WITH THE NORMALLY-GENERATED RETURNS !
%
%% Calculate Moments of the distributions :
% Empirical returns :
meanRet = mean(Ret);
varRet = var(Ret);
skewRet = skewness(Ret);
KurtRet = kurtosis(Ret);
% Normally-generated returns :
meanNorm = mean(NormDistrVar);
varNorm = var(NormDistrVar);
skewNorm = skewness(NormDistrVar);
KurtNorm = kurtosis(NormDistrVar);
%
% For example, we notice completely different KURTOSIS, different SKEWNESS,
% even though the mean and variance are the same.
%
%% Kolmogorov-Smirnov test :
% H0 = values belongs to normal distribution
x = (Ret-meanRet)/sqrt(varRet);
[h,pval] = kstest(x); % if h = 1 rejects the null hypothesis
%
% h = 1 : therefore we reject the null hypothesis : the empirical returns
% are not normally-distributed.
%
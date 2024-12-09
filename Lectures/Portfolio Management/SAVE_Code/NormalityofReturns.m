%% Are returns normal? 
clear all
close all
clc

%% Read Prices
path_map        = 'C:\Users\ginevra.angelini\Desktop\lezioni_poli\lezioni\Lezione3\';
filename        = 'geo_index_prices.xlsx';

table_prices = readtable(strcat(path_map, filename));
%% Calculate Returns (all period)
dt = table_prices(:,1).Variables;
values = table_prices(:,2).Variables;

Ret = values(2:end)./values(1:end-1);
% Plot
h = figure;
histogram(Ret, 'EdgeAlpha', 0.8)
%% Generation of N  random values from a normal distribution with same mean and variance
N =length(Ret);
NormDistrVar = zeros(1,N);
for i = 1:N
    r_ = random('Normal', mean(Ret), std(Ret));
    NormDistrVar(i) = r_;
end
%% Ensemble plot
h = figure;
histogram(Ret, 'EdgeAlpha', 0.8)
hold on 
histogram(NormDistrVar, 'EdgeAlpha', 0.7)
legend('Original Returtns Distr', 'Generated Values')

%% Calculate Moments of the distributions
meanRet = mean(Ret);
varRet = var(Ret);
skewRet = skewness(Ret);
KurtRet = kurtosis(Ret);

meanNorm = mean(NormDistrVar);
varNorm = var(NormDistrVar);
skewNorm = skewness(NormDistrVar);
KurtNorm = kurtosis(NormDistrVar);
%% Kolgornov-Smirnov H0 = values belongs to normal distribution

x = (Ret-meanRet)/sqrt(varRet);
[h,pval] = kstest(x); %if h = 1 rejects the null hypothesis

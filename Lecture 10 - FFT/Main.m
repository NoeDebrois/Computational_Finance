clear; close all;
param.rf=0.02; %risk-free rate
param.q= 0; %dividend
param.distr=1; %normal distribution
param.m=0; %drift
param.s=0.2; %volatility
S_0=1; K=1; Ndate=12; Barrier=0.8; N=2^14;
param.T=1; %maturity
param.dt=param.T/Ndate;
[S,v] = CONV( S_0, K, Ndate, N, Barrier, param);
price=interp1(S,v,S_0,'spline')
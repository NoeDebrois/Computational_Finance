clear
%% Pricing an European Call under VG
% Market/Contract Input
S0=100;
r=2/100;
K=105; T=2;
% Model (VG) Input
theta=0.2; sigma=0.6; k=0.5;
% MC Input
Nsim=1e6;

%% 1. Simulation
char_exp=@(u) -log(1+u.^2*sigma^2*k/2-1i*theta*k*u)/k;
drift=r-char_exp(-1i);
dS=k*icdf('gamma',rand(Nsim,1),T/k,1);
XT=0+drift*T+theta*dS+sigma*sqrt(dS).*randn(Nsim,1);
ST=S0*exp(XT);
disp('r-n check: must be zero')
[check,~,CI_check]=normfit(ST-S0*exp(r*T))
%% 2. Discounted Payoff
discpayoff=exp(-r*T)*max( ST-K,0);
%% 3. Price and Confidence Interval
[price,~,CI]=normfit(discpayoff)
priceCM=FFT_CM_Call_VG(K,[sigma,theta,k],T,r,S0)






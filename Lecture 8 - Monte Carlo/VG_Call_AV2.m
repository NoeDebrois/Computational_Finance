clear
%% Pricing an European Call under VG
%--> AV variance reduction
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
u=rand(Nsim,1);
dS1=k*icdf('gamma',u,T/k,1);
dS2=k*icdf('gamma',1-u,T/k,1);
z=randn(Nsim,1);
XT=0+drift*T+theta*dS1+sigma*sqrt(dS1).*z;
ST1=S0*exp(XT);
XT=0+drift*T+theta*dS2+sigma*sqrt(dS2).*z;
ST2=S0*exp(XT);
disp('r-n check: must be zero')
[check1,~,CI_check]=normfit(ST1-S0*exp(r*T))
[check2,~,CI_check]=normfit(ST2-S0*exp(r*T))
%% 2. Discounted Payoff
discpayoff1=exp(-r*T)*max( ST1-K,0);
discpayoff2=exp(-r*T)*max( ST2-K,0);
%% 3. Price and Confidence Interval
[price,~,CI]=normfit( (discpayoff1+discpayoff2)/2 )
priceCM=FFT_CM_Call_VG(K,[sigma,theta,k],T,r,S0)






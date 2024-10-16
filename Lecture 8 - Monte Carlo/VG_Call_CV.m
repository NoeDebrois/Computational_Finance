clear
%% Pricing an European Call under VG
%> CV variance reduction
% Market/Contract Input
S0=100;
r=2/100;
K=105; T=2;
% Model (VG) Input
theta=0.2; sigma=0.6; k=0.5;
% MC Input
Nsim=1e6; Nsim2=1000;

%% 1. Simulation
% f=ST
Ef=S0*exp(r*T);
char_exp=@(u) -log(1+u.^2*sigma^2*k/2-1i*theta*k*u)/k;
drift=r-char_exp(-1i);
%1. Small simulation to compute alpha
dS=k*icdf('gamma',rand(Nsim2,1),T/k,1);
XT=0+drift*T+theta*dS+sigma*sqrt(dS).*randn(Nsim2,1);
f=S0*exp(XT); g=exp(-r*T)*max( f-K,0);
A=cov(g,f);
alpha=-A(1,2)/A(2,2);
%2. Compute the price
dS=k*icdf('gamma',rand(Nsim,1),T/k,1);
XT=0+drift*T+theta*dS+sigma*sqrt(dS).*randn(Nsim,1);
f=S0*exp(XT); g=exp(-r*T)*max( f-K,0);
%% Price and Confidence Interval
[price,~,CI]=normfit(g+alpha*(f-Ef))
priceCM=FFT_CM_Call_VG(K,[sigma,theta,k],T,r,S0)






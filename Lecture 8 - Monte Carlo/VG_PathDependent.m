clear
%% Pricing Asian/Lookback under VG
% Market/Contract Input
S0=100;
r=2/100;
K=105; T=2;
% monitoring -> monthly
M=round(12*T);
% Model (VG) Input
theta=0.2; sigma=0.6; k=0.5;
% MC Input
Nsim=1e6;

%% 1. Simulation
dt=T/M;
X=zeros(Nsim,M+1);
char_exp=@(u) -log(1+u.^2*sigma^2*k/2-1i*theta*k*u)/k;
drift=r-char_exp(-1i);
dS=k*icdf('gamma',rand(Nsim,M),dt/k,1);
Z=randn(Nsim,M);
for j=1:M
    X(:,j+1)=X(:,j)+drift*dt+...
        theta*dS(:,j)+sigma*sqrt(dS(:,j)).*Z(:,j);
end
S=S0*exp(X); clear X dS Z;
disp('r-n check: must be zero')
[check,~,CI_check]=normfit(S(:,end)-S0*exp(r*T))
%% 2. Discounted Payoff
discpayoffEU=exp(-r*T)*max( S(:,end)-K,0);
% Asian Fixed Strike Put
discpayoffAFiP=exp(-r*T)*max( K-mean(S,2),0);
% Asian Floating Strike Call
discpayoffAFlC=exp(-r*T)*max( S(:,end)-mean(S,2),0);
% Lookback Fixed Strike Put on the minimum
discpayoffLFiPmin=exp(-r*T)*max( K-min(S,[],2),0);

%% 3. Price and Confidence Interval
[priceEU,~,CIEU]=normfit(discpayoffEU)
priceCMEU=FFT_CM_Call_VG(K,[sigma,theta,k],T,r,S0)
[priceAFiP,~,CIAFiP]=normfit(discpayoffAFiP)
[priceAFlC,~,CIAFlC]=normfit(discpayoffAFlC)
[priceLFiPmin,~,CILFiPmin]=normfit(discpayoffLFiPmin)





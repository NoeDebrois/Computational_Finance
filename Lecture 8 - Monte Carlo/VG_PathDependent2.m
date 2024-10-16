clear
%% Pricing Knock&Out Call Option under VG
% Market/Contract Input
S0=100;
r=2/100;
K=105; T=2; U=130; L=80;
% monitoring -> monthly
M=round(12*T);
% Model (VG) Input
theta=0.2; sigma=0.6; k=0.5;
% MC Input
Nsim=1e5;

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
discpayoff=exp(-r*T)*max( S(:,end)-K,0)...
    .*(min(S,[],2)>L).*(max(S,[],2)<U);
%% 3. Price and Confidence Interval
[price,~,CI]=normfit(discpayoff)





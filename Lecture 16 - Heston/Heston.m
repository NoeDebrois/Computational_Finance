% Simulate the Heston model via the Euler scheme
clear; close all;
% Input
%>variance
V0=0.1; % variance at time 0
theta=0.3; % mean 
k=0.1; % speed of mean reversion
epsilon=0.2; %vol-of-var
%>price
r=0.02;
S0=100; 
%>correlation between price and variance
rho=-0.6;
% Calibration
%> r, S0 observed on the market
%> [V0, theta, k, epsilon, rho] have to be calibrated (from European Call
%exploiting Carr-Madan --> characteristic function of Heston known)

Feller_condition=( epsilon^2-2*k*theta )
if Feller_condition<=0
    disp('Feller Condition is satisfied')
else
    disp('Feller Condition is NOT satisfied')
end
Nsim=1e4; T=1; Nsteps=25; %discretization
% MC simulation
dt=T/Nsteps; t=linspace(0,T,Nsteps+1); % time grid
x=zeros(Nsim,Nsteps+1); x(:,1)=log(S0); %log price process
V=zeros(Nsim,Nsteps+1); V(:,1)=V0; %variance process
% mean and variance-covariance matrix of the normal distribution to be sampled
mu=[0;0]; VC=[1 rho;rho 1]; 
for i=1:Nsteps
    Z=mvnrnd(mu,VC,Nsim); 
    x(:,i+1)=x(:,i)+(r-V(:,i)/2)*dt+sqrt(V(:,i)*dt).*Z(:,1);
    V(:,i+1)=V(:,i)+k*(theta-V(:,i))*dt+...
        epsilon*sqrt(V(:,i)*dt).*Z(:,2);
end
S=exp(x);
figure; hold on
subplot(1,2,1); plot(t,S); title('Price');
subplot(1,2,2); plot(t,V); title('Variance');


    
    




clear; close all;
%% Simulation of the Merton model
mu=0.05; %drift
sigma=0.4; %volatility
lambda=2; %Poisson intensity
muJ=0.01; deltaJ=0.2; %Jumpsize parameters

T=2 % Maturity
M=100; % Number of steps in time
dt=T/M;
X=zeros(M+1,1);
Z=randn(M,1);
for i=1:M
    X(i+1)=X(i)+mu*dt+sigma*sqrt(dt)*Z(i);
    % jumps term
    Ndt=poissrnd(lambda*dt);
    if Ndt==0
        J=0;
    else
        J=sum(muJ+deltaJ*randn(Ndt,1));
    end
    X(i+1)=X(i+1)+J;
end
plot(exp(X))

    
    
    
    
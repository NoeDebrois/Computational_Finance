clear
close all
% parameters
y0=0.02;
beta=0.3; lambda=1; eta=0.05;
if beta^2-2*lambda*eta<=0
    disp('Feller condition satisfied')
else
    disp('Feller condition NOT satisfied')
end
T=1; M=100; Nsim=100; 
dt=T/M;
y=zeros(Nsim,M+1); y(:,1)=y0;
for i=1:M
    y(:,i+1)=y(:,i)+lambda*(eta-y(:,i))*dt+...
             beta*sqrt(y(:,i))*sqrt(dt).*randn(Nsim,1);
    y(:,i+1)=max(y(:,i+1),0);     
end
t=linspace(0,T,M+1); plot(t,y);






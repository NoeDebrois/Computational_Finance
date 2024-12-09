B=0.98;
t=1;
f=1;
K=1.1;
%% HJM modelling of swaps p=2,n=0, price european call
T_1=2;
a=0.1;
sigma_1= 0.2; %*exp(a*(T_1-t));
sigma_2=0.1;

d2=(log(f/K)-0.5*(sigma_2^2*t+sigma_1^2/a/4*(exp(2*a*T_1)-exp(2*a*(T_1-t)))))/sqrt(sigma_2^2*t+sigma_1^2/a/4*(exp(2*a*T_1)-exp(2*a*(T_1-t))));
d1=d2+sqrt(sigma_2^2*t+sigma_1^2/a/4*(exp(2*a*T_1)-exp(2*a*(T_1-t))));
B*(f*normcdf(d1)-K*normcdf(d2))
%% HJM modelling p=0,n=1 (Driven by Levy NIG)
eta=0.2;
Strike=1.01
% parameter
F0=1; T=1; B=0.96; sigma=0.17801;
k_IG=0.5; theta=-0.5;
% discretization parameter
Npow=20;  A=1200;
V=@(v) T*(1/k_IG)*(1-sqrt(1+v.^2*sigma^2*k_IG*eta^2-2*1i*theta*v*k_IG*eta));
CharFunc=@(v) exp(V(v)-1i*v*V(-1i));
FFT_CM_Call(Strike,F0,B,CharFunc,Npow,A)
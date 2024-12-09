clear all

%% Simulation: Arithmetic model framework p=1, m=1, n=1;
%Driven by Gamma process.
t=1;
dt=t/200;
time_steps=t/dt;
%% Seasonality
Lambda=@(t) 1+sin(4*pi*t);
%% Constant Parameters
a=0.5;
mu=0;
sigma=0.2;
b=1;
d=0;
eta=0.1;
x_0=1;
y_0=1;
%% Gamma parameter
k=0.5;
shape=1/k*dt; scale=k;
nu=@(y) (y>1e-8)*1/k.*exp(-y/k)./y;
gamma=quadgk(@(y) nu(y).*y,-1,1);
%% Simulate OU spot
N_sim=10^5;
path = gamrnd(shape,scale,N_sim,time_steps).*exp(-b*((t:-dt:dt)-dt/2)); 

G=randn(N_sim,1);
X=x_0*exp(-a*t)+mu/a*(1-exp(-a*t))+sigma*sqrt((1-exp(-2*a*t))/2/a)*G;
Y=y_0*exp(-b*t)+d/b*(1-exp(-b*t))+eta*sum(path,2);
S= Lambda(t)+X+Y;
mean(S)
% Don't forget to check the condition A!
%% Simulate fwd under Q
theta_x=@(t) 0;
theta_y=@(t) 0;

t=1;
T=2;

Theta=@(t,T)  eta*quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y)*exp(-b*(T-u)).*y.*(exp(theta_y(u)*y)-(y<1)),1e-8,30),u),t,T)...
    +gamma/b*eta*(1-exp(-b*(T-t)))+sigma*quadgk(@(u)theta_x(u)*exp(-a*(T-u)),t,T);
X_Q=X;
Y_Q=Y-gamma/b*eta*(1-exp(-b*t))-eta*quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y).*y*exp(-b*u),1,30),u),0,t);
f=@(t,T,X,Y) Lambda(T)+Theta(t,T)+mu/a*(1-exp(-a*(T-t)))+d/b*(1-exp(-b*(T-t)))+X_Q*exp(-a*(T-t))+Y_Q*exp(-b*(T-t));
mean(f(t,T,X_Q,Y_Q))
%% Simulate swap under Q (Proposition 4.14)

T1=2;
T2=3;
%We have the closed formula!
F=@(t,T1,T2,X,Y) quadgk(@(u) arrayfun(@(u)Theta(t,u),u),T1,T2)+quadgk(@(t)Lambda(t),T1,T2)/(T2-T1)+X_Q*(exp(-a*T1)-exp(-a*T2))/(a*(T2-T1))+Y_Q*(exp(-b*T1)-exp(-b*T2))/(b*(T2-T1));
mean(F(t,T1,T2,X_Q,Y_Q))

%% Pricing option 

mean(max(F(t,T1,T2,X,Y)-1.2,0))


 

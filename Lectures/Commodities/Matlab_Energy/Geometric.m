%% Simulation: geometric model framework p=1, m=1, n=1;
%Driven by Gamma process.
t=1;
dt=t/200;
time_steps=t/dt;
%% Seasonality (this is a simple one)
Lambda=@(t) 1+sin(2*pi*t);
%% Constant Parameters
a=0.5;
mu=0;
sigma=0.2;
b=1;
d=0;
eta=0.1;
%starting points
x_0=1;
y_0=1;
%% Gamma parameter
k=0.5;
shape=1/k*dt; scale=k; %Writing of the gamma parameters with respect to a single quantity k
nu=@(y) (y>1e-8)*1/k.*exp(-y/k)./y;
gamma=quadgk(@(y) nu(y).*y,-1,1);
%Check condition GG!!
%% Simulate Simple Levy case
N_sim=10^5;

gg= @(z) exp(t*(1i*z*gamma+quadgk(@(y) nu(y).*(exp(z*1i*y)-1-(z*1i*y).*(abs(y)<1)),1e-8,50))); %char function

path=gamrnd(shape,scale,N_sim,time_steps); %simulating increments

YY=sum(path,2);

mean(YY)
%check of martingality
mean(exp(YY))
gg(-1i)

%% Simulate OU (Spot)
%path contains inrements
path = gamrnd(shape,scale,N_sim,time_steps).*exp(-b*((t:-dt:dt)-dt/2)); % To check dt

G=randn(N_sim,1);
X=x_0*exp(-a*t)+mu/a*(1-exp(-a*t))+sigma*sqrt((1-exp(-2*a*t))/2/a)*G;
Y=y_0*exp(-b*t)+d/b*(1-exp(-b*t))+eta*sum(path,2); %With sum, vector N_sim x 1
S=Lambda(t)*exp(X+Y);

phi_x=@(z) exp(1i*z*(x_0*exp(-a*t)+(mu/a)*(1-exp(-a*t)))-1/(4*a)*z^2*sigma^2*(1-exp(-2*a*t)));
phi_y=@(z) exp(1i*z*(y_0*exp(-b*t)+(d/b)*(1-exp(-b*t)))+1i*z*gamma/b*eta*(1-exp(-b*t))+...
    quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y).*(exp(z*1i*eta*y*exp(-b*(t-u)))-1-z*1i*eta*y.*exp(-b*(t-u)).*(abs(y)<1)),1e-8,30),u),0,t));
%Check
phi_x(-1i)
mean(exp(X))

phi_y(-1i)
mean(exp(Y))

%% Simulate Fwd under Q
%Let suppose parameters of esscher transform equal to zero
theta_x=@(t) 0;
theta_y=@(t) 0;

t=1;
T=2;

Theta=@(t,T) exp(  quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y).*(exp(eta*y*exp(-b*(T-u))+theta_y(u))-1-(eta.*exp(-b*(T-u))+theta_y(u))*y.*(abs(y)<1)),1e-8,30),u),t,T)...
    -quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y).*(exp(theta_y(u)*y)-1-theta_y(u)*y*exp(-b*(T-u)).*(abs(y)<1)),1e-8,30),u),t,T)...
    +gamma/b*eta*(1-exp(-b*(T-t))+1/4*sigma^2/a*(1-exp(-2*a*(T-t)))+sigma*quadgk(@(u)theta_x(u)*exp(-a*(T-u)),t,T)));
%Change in the drift for esscher
phi_XX=@(z) exp(-1/(4*a)*z^2*sigma^2*(1-exp(-2*a*t)));
phi_YY=@(z) exp(1i*z*gamma/b*eta*(1-exp(-b*t))+quadgk(@(u)arrayfun(@(u)quadgk(@(y) nu(y).*(exp(z*1i*eta*y*exp(-b*(t-u)))-1-z*1i*eta*y.*exp(-b*(t-u)).*(abs(y)<1)),1e-8,30),u),0,t));
X_Q=X*phi_XX(-1i);
Y_Q=Y*phi_YY(-1i);
f=@(t,T,X_Q,Y_Q) Lambda(T)*Theta(t,T).*exp(mu/a*(1-exp(-a*(T-t)))+d/b*(1-exp(-b*(T-t)))+exp(-a*(T-t))*X_Q+exp(-b*(T-t))*Y_Q);
%The estimated price of the forward by simulation
mean(f(t,T,X_Q,Y_Q))
%% Simulate Swap under Q (settlement at the end of delivery period).

T1=2;
T2=3;
Times=T1:((T2-T1)/12):T2; %Suppose monthly approx.
F=0;
for i=2:length(Times)
F=F+f(t,Times(i),X_Q,Y_Q)/(T2-T1)*(Times(i)-Times(i-1));
end
%Why we need this approx? It is a geometric model!
mean(F)


%% Price Option on Swap under Q (maturity at time t).
B=0.98; %discount factor
K=1.8;
C=B*mean(max(F-K,0));
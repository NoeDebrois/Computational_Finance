clear
close all
% Pricing a plain vanilla Call Option
% B&S model - logprice PDE
% Implicit Euler

%% 1. Parameters & Grids
S0=1; K=0.95; r=0.001; sigma=0.5; T=1;
M=1000;
N=2000;
Smin=0.1*S0; Smax=3*S0;
xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
%x=linspace(xmin,xmax,N+1);
x=xmin+(0:N)*dx;
dt=T/M;
t=(0:M)*dt;

%% Implicit Euler - Matrix Construction
A=dt*(-(r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
B=-1+dt*(-sigma^2/(dx^2)-r);
C=dt*( (r-sigma^2/2)/(2*dx)+sigma^2/(2*dx^2));
Mat=spalloc(N+1,N+1,3*(N-1)+2);
Mat(1,1)=1;
for i=2:N
    Mat(i, [i-1,i,i+1])=[A B C];
end
Mat(end,end)=1;

%% Implicit Euler - Backward in time loop
c=max(S0*exp(x')-K,0); 
rhs=zeros(N+1,1);
for j=M:-1:1
    % know c t_j -> compute cnew t_{j-1}
    rhs(1)=0;
    rhs(2:end-1)=-c(2:end-1);
    rhs(end)=Smax-K*exp(-r*(T-(j-1)*dt));
    c=Mat\rhs;
end
figure
plot(S0*exp(x),c); title('Call price'); xlabel('S');
price=interp1(x,c,0,'spline')
price_ex=blsprice(S0,K,r,T,sigma)





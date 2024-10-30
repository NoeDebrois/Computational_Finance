clear
close all
% Pricing a plain vanilla Call Option
% B&S model - logprice PDE
% Explicit Euler

%% 1. Parameters & Grids
S0=1; K=0.95; r=0.001; sigma=0.5; T=1;
M=1000;
N=300;
Smin=0.1*S0; Smax=3*S0;
xmin=log(Smin/S0); xmax=log(Smax/S0);
dx=(xmax-xmin)/N;
%x=linspace(xmin,xmax,N+1);
x=xmin+(0:N)*dx;
dt=T/M;
t=(0:M)*dt;

%% Explicit Euler - Backward in time loop
c=max(S0*exp(x')-K,0);
cnew=zeros(size(c));
for j=M:-1:1
    % know c t_j -> compute cnew t_{j-1}
    cnew(1)=0;
    for i=2:N
        cnew(i)=c(i)+dt*( ...
            (r-sigma^2/2)*(c(i+1)-c(i-1))/(2*dx)+...
            sigma^2/2*(c(i+1)-2*c(i)+c(i-1))/dx^2+...
            -r*c(i) );
    end
    cnew(N+1)=Smax-K*exp(-r*(T-(j-1)*dt));
    c=cnew;
end
figure
plot(S0*exp(x),c); title('Call price'); xlabel('S');
price=interp1(x,c,0,'spline')
price_ex=blsprice(S0,K,r,T,sigma)





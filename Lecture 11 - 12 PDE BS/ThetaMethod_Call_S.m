clear
close all
% Pricing a plain vanilla Call Option
% B&S model - price PDE
% ThetaMethod

%% 1. Parameters & Grids
S0=1; K=0.95; r=0.001; sigma=0.5; T=1;
M=1000;
N=2000;
theta=0.5; % 0 Implicit, 0.5 CN, 1 Explicit
Smin=0.1*S0; Smax=3*S0;
dS=(Smax-Smin)/N;
S=Smin+(0:N)*dS;
dt=T/M;
t=(0:M)*dt;

%% Matrix Construction
A=@(S)(1-theta)*dt*(-r*S/(2*dS)+sigma^2*S^2/(2*dS^2));
B=@(S)-1+dt*(1-theta)*(-sigma^2*S^2/(dS^2)-r);
C=@(S)dt*(1-theta)*( r*S/(2*dS)+sigma^2*S^2/(2*dS^2));
Mat=spalloc(N+1,N+1,3*(N-1)+2);
Mat(1,1)=1;
for i=2:N
    Mat(i, [i-1,i,i+1])=[A(S(i)) B(S(i)) C(S(i))];
end
Mat(end,end)=1;
Ah=@(S)-(theta)*dt*(-r*S/(2*dS)+sigma^2*S^2/(2*dS^2));
Bh=@(S)-1-dt*(theta)*(-sigma^2*S^2/(dS^2)-r);
Ch=@(S)-dt*(theta)*( r*S/(2*dS)+sigma^2*S^2/(2*dS^2));
Mat_rhs=spalloc(N+1,N+1,3*(N-1));
for i=2:N
    Mat_rhs(i, [i-1,i,i+1])=[Ah(S(i)) Bh(S(i)) Ch(S(i))];
end

%% Backward in time loop
c=max(S'-K,0); 
rhs=zeros(N+1,1);
for j=M:-1:1
    % know c t_j -> compute cnew t_{j-1}
    rhs=Mat_rhs*c;
    rhs(1)=0;
    rhs(end)=Smax-K*exp(-r*(T-(j-1)*dt));
    c=Mat\rhs;
end
figure
plot(S,c); title('Call price'); xlabel('S');
price=interp1(S,c,S0,'spline')
price_ex=blsprice(S0,K,r,T,sigma)





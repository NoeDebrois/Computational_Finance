function [Error,Price]=fun_CM_BS(x,spot,strike,rf,maturity,pmkt)

% Initialize
Price=zeros(length(strike),1);

% Delete duplicate
mat=unique(maturity);

for i=1:length(mat)
    index=find(maturity==mat(i));
    Price(index)=FFT_CM_Call(strike(index),spot,mat(i),rf,x);
end
Error=abs(pmkt-Price);

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

function [Price]=FFT_CM_Call(Strike,S0,T,r,sigma)
% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: BS

% discretization parameter
Npow=15;
N=2^Npow;
A=600;
% v-> compute integral as a summation
eta=A/N;
v=[0:eta:A*(N-1)/N]; v(1)=1e-22;
% lambda-> compute summation via FFT
lambda=2*pi/(N*eta);
k=-lambda*N/2+lambda*(0:N-1);

% Fourier transform of z_k
CharFunc=@(v) exp(T*CharExpBS(v,sigma));
Z_k=exp(1i*r*v*T).*(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));
% Option Price
w=ones(1,N); w(1)=0.5; w(end)=0.5;
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);
C=S0*(z_k+max(1-exp(k-r*T),0));
K=S0*exp(k);

% Output
index=find( K>0.1*S0 & K<3*S0 );
C=C(index); K=K(index);
Price=interp1(K,C,Strike,'spline');

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
function V=CharExpBS(v,sigma)
% risk-neutral characteristic exponent
V=-sigma^2/2*1i*v-sigma^2/2*v.^2;

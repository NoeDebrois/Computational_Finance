function [Price]=FFT_CM_Call(Strike,F0,B,CharFunc,Npow,A)
% Price of a Plain Vanilla Call exploiting the Carr-Madan algorithm
% Model: BS
format long


% discretization parameter
 N=2^Npow;

% v-> compute integral as a summation
eta=A/N; v=[0:eta:A*(N-1)/N]; v(1)=1e-22;
% lambda-> compute summation via FFT
lambda=2*pi/(N*eta); 
k=-lambda*N/2+lambda*(0:N-1);

% Fourier transform of z_k

%CharFunc=@(v) exp(T*CharExpBS(v,sigma));

disp('RiskNeutral Check')
CharFunc(-1i)
Z_k=(CharFunc(v-1i)-1)./(1i*v.*(1i*v+1));
% Option Price
w=ones(1,N); w(1)=0.5; w(end)=0.5;
x=w.*eta.*Z_k.*exp(1i*pi*(0:N-1));
z_k=real(fft(x)/pi);
C=B*F0*(z_k+max(1-exp(k),0));
K=F0*exp(k);

% Output
index=find( K>0.1*F0 & K<3*F0 );
C=C(index); K=K(index);
plot(K,C)
title( 'Option Price' );
xlabel('Strike');
Price=interp1(K,C,Strike,'spline');
disp('Check')

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>










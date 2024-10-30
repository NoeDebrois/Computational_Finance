function [S,v] = CONV( S_0, K, Ndate, N,...
    Barrier, param)

b=2.5;
[x,~,~,H] = kernel(N,-b,b,param,0);
S = S_0*exp(x);
% Payoff and Fourier transform of the density
v = max(S-K,0).*(S>Barrier);
H=ifftshift(H);
for j = Ndate:-1:1
    v=real((fft(ifft(v).*H)))*exp(-param.rf*param.dt);
%     v=real(fftshift(fft...
%         (ifft(ifftshift(v)).*H)))*...
%         exp(-param.rf*param.dt);
    v(S<=Barrier) = 0; 
end
index=find( (S>0.1*S_0).*(S<3*S_0));
S=S(index); v=v(index);
figure
plot(S,v);
xlabel('S'); title('Call Option');


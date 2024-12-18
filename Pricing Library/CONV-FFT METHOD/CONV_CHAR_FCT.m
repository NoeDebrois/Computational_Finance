function [S,v] = CONV_CHAR_FCT(S_0, K, Ndate, N, Barrier, param)

b = 2.5;
[x, ~, ~, H] = KERNEL_CHAR_FCT(N, -b, b, param, 0);
S = S_0 * exp(x);

v = max(S - K, 0) .* (S > Barrier);
H = ifftshift(H);

for j = Ndate:-1:1
    v = real((fft(ifft(v) .* H))) * exp(-param.rf * param.dt);
    v(S <= Barrier) = 0;
end

index = find((S > 0.1 * S_0) .* (S < 3 * S_0));
S = S(index); v = v(index);
figure
plot(S, v);
xlabel('Grid of spot price S'); 
ylabel('Call Option Price (v)')
title('Barrier (Down&Out) Call Option Price (v) vs a grid of Spot Price (S)');

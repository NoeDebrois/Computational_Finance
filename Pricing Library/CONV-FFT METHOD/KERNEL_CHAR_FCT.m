function [x, h, w, H] = KERNEL_CHAR_FCT(ngrid, xmin, xmax, parameters, flag)

if nargin == 4
    flag = 0;
end

N = ngrid / 2;
dx = (xmax - xmin) / ngrid;
x = dx * (-N:N-1);

dw = 2 * pi / (xmax - xmin);
w = dw * (-N:N-1);

H = CHAR_FCT(w, parameters, flag);

h = real(fftshift(fft(ifftshift(H)))) / (xmax - xmin);

figure
plot(x, h)
title("Density of -X, on the log-price grid (x)");
xlabel("Log-price grid (x)");
ylabel("Density of -X")

figure
plot(w, abs(H))
title("Conjugated characteristic function of X, on the Fourier space grid (w)");
xlabel("Fourier space grid (w)");
ylabel("Conjugated characteristic function of X")

end
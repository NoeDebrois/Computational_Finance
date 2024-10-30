function [x,h,w,H] = kernel(ngrid,xmin,xmax,parameters,flag)
% flag=0 --> characteristic function for backward-in-time 
% flag=1 --> characteristic function for forward-in-time

if nargin==4
    flag=0; %backward characteristic function
end
    
N = ngrid/2;
dx = (xmax-xmin)/ngrid;
x = dx*(-N:N-1);
dw = 2*pi/(xmax-xmin);
w = dw*(-N:N-1);

H = charfunction(w,parameters,flag); % characteristic function
h = real(fftshift(fft(ifftshift(H))))/(xmax-xmin); % kernel

figure
plot(x,h)
figure
plot(w,H)
% Basically, MATLAB implements an fft (as do other sets of functions) 
% whose results AND inputs are ordered (for an array of N elements), 
% from 0 to (N/2-1) and then from –N/2 to -1. 
% (I will call this 'swapped' from now on.) 
% But in what I might call a 'natural' ordering, 
% for a spectrum or time function centred at zero and extending equal time 
% (or frequency) either side of zero, 
% the arrays of input data and results would be ordered (2’s complement style)
% from –N/2 to N/2-1 where N is the number of elements. 
% That is, time zero (or frequency zero) are at the centre of such an array. 
% 
% It is well known that the results of MATLAB's fft() function must be 
% re-ordered with the MATLAB function fftshift() so that they are 'naturally' 
% ordered (for instance when plotting).
% 
%    H = fftshift( fft ( h ) );
% 
 

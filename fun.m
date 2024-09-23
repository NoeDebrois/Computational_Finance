function dp=fun(vol,spot,strike,rf,maturity,pmkt)
% Compute B&S Price according to the B&S formula :
% d1 :
d1 = (log(spot./strike)+(rf+vol.*vol/2).*maturity)./(vol*sqrt(maturity));
% d2 :
d2 = d1 - vol.*sqrt(maturity);
% N(d1) :
nd1 = normcdf(d1,0,1);
% N(d2) :
nd2 = normcdf(d2,0,1);
% B&S price according to the B&S formula :
bs = (spot.*nd1 - strike.*exp(-rf.*maturity).*nd2);

% Compute the difference between the B&S price and the real price :
dp = bs - pmkt; % This is what we want to minimize when calibrating.
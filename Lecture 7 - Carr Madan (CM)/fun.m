function dp=fun(params,spot,strike,rf,maturity,pmkt)
T=unique(maturity); Price=ones(size(strike));
for i=1:length(T)
    idx=find(maturity==T(i));
    Price(idx)=FFT_CM_Call_Kou...
        (strike(idx),params,T(i),rf,spot);
end
%compute the difference between the
%model price and the real price
dp=Price-pmkt;
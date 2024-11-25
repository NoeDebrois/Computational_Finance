function dp = fun_CM_BS(params, spot, strike, rf, maturity, pmkt)
% Computes the difference between the model (here : B&S) price and the real
% price.
%
% CM algorithm works for different strikes on the same maturity :
T = unique(maturity); % Set of maturities that we work with (no double).
Price = ones(size(strike)); % One price per strike.


for i=1:length(T)
    idx = find(maturity == T(i));

    Price(idx) = FFT_CM_Call_BS(strike(idx), params, T(i), rf, spot);
end

% Difference between model price and market (real) price :
dp = Price - pmkt;
%% LECTURE 7 - Carr-Madan Method - Noé Debrois - 27/10/2024
% This code implements the computation of the difference between the price
% under Carr-Madan (CM) algorithm using the Merton model for jump-diffusion
% processes (model price) and the price of the market (real price). This 
% function is used in order to calibrate our model to the market. This
% function call FFT_CM_Call_Merton.m function to compute the model price.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dp = fun_MERTON(params, spot, strike, rf, maturity, pmkt)
% Computes the difference between the model price and the real price.
%
% CM algorithm works for different strikes on the same maturity :
T = unique(maturity); % Set of maturities that we work with (no double).
Price = ones(size(strike)); % One price per strike.

% Compute the model price using CM :
for i=1:length(T)
    % For each maturity in our set of (unique) maturities, we compute the 
    % CM price using all the strikes. We will, at the end, have a vector 
    % with 'length(T)' elements (one per unique maturity). Each element
    % will be a vector of prices for THAT maturity (only one), and THOSE 
    % strikeS (if many).
    %
    idx = find(maturity == T(i)); % find() is locating the indices in the
    % maturity array where the elements are equal to T(i).
    Price(idx) = FFT_CM_Call_Merton(strike(idx), params, T(i), rf, spot);
    % We use the strikes corresponding to that maturity (CM : 1 maturity, 
    % different strikes).
end

% Difference between model price and market (real) price :
dp = Price - pmkt;
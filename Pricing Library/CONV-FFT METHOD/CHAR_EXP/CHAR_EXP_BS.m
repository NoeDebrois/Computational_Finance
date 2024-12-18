function V = CHAR_EXP_BS(v, parameters)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returns the characteristic exponent of logprices in BS model
    % INPUT: v = point of evaluation; 
    %        parameters = parameters struct{
    %                     - sigma = volatility of the model
    %                     }
    % OUTPUT: characteristic exponent for the BS model under risk neutral measure 
    % cf the fact that for N(mu, sigma^2) : CF(u) = exp(i*mu*u - sigma^2/2 * u^2).
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % De-struct parameters
    sigma = parameters.sigma;
    
    % Evaluate char exp under risk-neutral measure 
    V = @(v) - sigma^2/2 * v.^2; % Without drift
    drift_rn = -V(-1i); % Drift Risk_neutral
    V = drift_rn * 1i * v + V(v);
end
%% LECTURE 10 - Convolution / FFT Methods - Noé Debrois - 30/10/2024
% This code implements the computation of characteristic function for a
% given distribution (case 1 : Normal, case 2 : Normal Inverse Gaussian).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = charfunction(u, parameters, flag)
% Computes the characteristic function (conjugated OR not) of a given
% distribution (Normal or Normal Inverse Gaussian) with a correction to
% take into account the RISK-FREE RATE, the potential DIVIDENDS, and the
% fact that, under Q measure, phi_X(-i) = 1 (i.e Psi_X(-i) = 0).
%
% cf FFT.pdf & CONV.pdf : in the computations, to make a convolution appear
% we need to complex-conjugate the characteristic function. That's why :
% flag = 0 for complex conjugate -> characteristic fct backward-in-time ; 
% flag = 1 for NO complex conjugate -> characteristic fct forward-in-time.

if nargin == 2
    flag = 0; % By default, we need the complex conjugate.
end

% To make it risk-neutral :
meancorrection = (parameters.rf - parameters.q) * parameters.dt...
-log(charfunction0(-1i, parameters));

F = exp(1i * meancorrection * u) .* charfunction0(u, parameters);

if flag == 0
    F = conj(F); % We need the complex conjugate to have a convolution.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = charfunction0(u, parameters)
% Computes the characteristic function of a Normal distribution (case 1) or
% a Normal Inverse Gaussian distribution (case 2).

dt = parameters.dt;
switch parameters.distr

    case 1 % Normal distribution.
        %
        % CF(u) = exp(i * mu * u - sigma^2 * t^2 / 2).
        %

        % Parameters :
        m = parameters.m;
        s = parameters.s;
	
	    % Rearrange parameters (time rescaling) :
        m = m * dt;
        s = s * sqrt(dt);

        % Computation of the CF (cf formula above) :
        F = exp(1i * u * m - 0.5 * (s * u).^2);

    case 2 % Normal inverse Gaussian (NIG) distribution.
        %
        % CF =
        %

        % Parameters :
        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta;

        % Rearrange parameters (time rescaling)
        alpha = alpha;
        beta = beta;
        delta = delta * dt;

        % Computation of the CF (cf formula above) :
        F = exp(-delta * (sqrt(alpha^2 - (beta + 1i * u).^2)...
            - sqrt(alpha^2 - beta^2)));
   
end

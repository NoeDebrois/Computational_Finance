function F = charfunction(u,parameters,flag)
% flag=0 --> characteristic function for backward-in-time 
% flag=1 --> characteristic function for forward-in-time

if nargin==2
    flag=0;
end

meancorrection = ...
(parameters.rf-parameters.q)*parameters.dt...
-log(charfunction0(-1i,parameters));
F = exp(1i*meancorrection*u).*...
       charfunction0(u,parameters);
if flag==0
    F=conj(F);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = charfunction0(u,parameters)
dt=parameters.dt;
switch parameters.distr

    case 1 % Normal
        m = parameters.m;
        s = parameters.s;
	
	    % Rearrange parameters (time rescaling)
        m = m*dt;
        s = s*sqrt(dt);

        F = exp(1i*u*m-0.5*(s*u).^2);

    case 2 % Normal inverse Gaussian (NIG) 
        alpha = parameters.alpha;
        beta = parameters.beta;
        delta = parameters.delta;

        % Rearrange parameters (time rescaling)
        alpha = alpha;
        beta = beta;
        delta = delta*dt;

        F = exp(-delta*(sqrt(alpha^2-(beta+1i*u).^2)-sqrt(alpha^2-beta^2)));
   
end

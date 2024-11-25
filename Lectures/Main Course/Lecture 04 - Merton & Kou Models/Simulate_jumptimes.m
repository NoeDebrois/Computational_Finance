% Two ways to simulate the Poisson Process and when the jumps occur.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%% Parameters :
%
lambda = 2; % for the Poisson distribution (don't forget to * by t).
t = 1;
% The rate parameter for the Poisson distribution @ t is lambda*t.
%
%% 1) Conditional Simulation : the easiest way.
%
% Generates random nb from Poisson distribution of rate parameter lambda*t:
Nt = poissrnd(lambda * t) %;
% Same as :
% Nt = icdf('Poisson',rand,lambda*t) where rand returns a RV ~ U([0, 1]).
%
% Jump times ordered chronologically :
T = sort(rand(1, Nt) * t) %;
%
%% 2) Countdown Simulation : a little bit harder.
%
T=[]; % Jump times
T(1) = 0; % We start @ 0
while T(end) < t % While we did not reach the "final time"...
    % Simulate/Sample the interarrival time ~ Exp(lambda) :
    tau = exprnd(1/lambda); 
    % Same as (cf help) :
    % icdf('Exponential', rand, 1/lambda) : where we must pass 1/lambda).
    %
    % Compute jump time :
    T = [T, T(end)+tau]; % Jump happens @ previous time + interarrival time
end
T = T(2:end-1) % The 1st elt is 0 but it's not a jump, & the last is after
% t, so we must remove it.
Nt = length(T) % Number of jumps.
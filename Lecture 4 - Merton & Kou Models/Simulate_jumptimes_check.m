% Two ways to simulate the Poisson Process and when the jumps occur.
% We simulate a lot of times with the 2 methods and compare mean & var.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
%% Parameters :
%
lambda = 2; % for the Poisson distribution (don't forget to * by t).
t = 1;
% The rate parameter for the Poisson distribution @ t is lambda*t.
Nsim = 1e6; % Number of simulations.
%
%% 1) Conditional Simulation : the easiest way.
% Create a vector of Poisson RV of length Nsim, with rate lambda*t :
Nt = poissrnd(lambda*t,Nsim,1); % Same as Nt=icdf('Poisson',rand,lambda*t);
% Print the mean and the variance of the number of jumps vector :
disp("Conditional simulation results : [mean, variance]")
mean(Nt)
var(Nt)
%
%% 2) Countdown Simulation : a little bit harder.
%
vN = zeros(Nsim,1);
%
for i=1:Nsim
    T = []; % Jump times
    T(1) = 0; % We start @ 0
    while T(end) < t % While we did not reach the "final time"...
        % Simulate/Sample the interarrival time :
        tau = exprnd(1/lambda); % Same as icdf('Exponential',rand,1/lambda) : 
        % check the help about how the pdf is implemented (we must pass
        % 1/lambda).
        %
        % Compute jump time :
        T = [T, T(end)+tau]; % Jump happens @ previous time + interarrival 
        % time
    end
    T = T(2:end-1); % The 1st elt is 0 but it's not a jump, & the last is
    % after t, so we must remove it.
    Nt = length(T); % Number of jumps.
    vN(i) = Nt; % Number of jumps for the i-th simulation.
end
% Print the mean and the variance of the number of jumps vector :
disp("Countdown simulation results : [mean, variance]")
mean(vN)
var(vN)
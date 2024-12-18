% Test parameters for Black-Scholes model
parameters.s = 0.2;  % Volatility
parameters.m = 0.0; 
parameters.dt = 1/252;   % Time step (assuming daily data with 252 trading days per year)
parameters.distr = 1;    % 1 corresponds to Normal distribution (Black-Scholes)
parameters.rf = 0.05;    % Risk-free rate
parameters.q = 0.02;     % Dividend yield

% Define the grid size and range
ngrid = 512;         % Number of grid points
xmin = -5;           % Minimum value of log-price
xmax = 5;            % Maximum value of log-price

% Call the kernel function for characteristic function and density
[x, h, w, H] = kernel(ngrid, xmin, xmax, parameters, 0);

% Display results
disp('Density of -X on the log-price grid (x):');
disp(h);

disp('Magnitude of the characteristic function of X on the Fourier space grid (w):');
disp(abs(H));

function [xbar, computing_time] = PrimalDualMin(x0, prox_fstar, prox_g, A, Astar, options)

if isempty(x0)
    disp("PrimalDualMinError:x0 is empty ... program about to crash!")
end

iter_max = 2000;
theta = 1; 
sigma = 0.5; 
tau = 0.01; 

if ~(nargin < 6 || isempty(options))
    if isfield(options, 'iter_max')
        iter_max = options.iter_max;
    end
    if isfield(options, 'theta')
        theta = options.theta;
    end
    if isfield(options, 'sigma')
        sigma = options.sigma;
    end
    if isfield(options, 'tau')
        tau = options.tau;
    end
end


[m,n] = size(x0);

% Some initialiization
x =x0;
xi = A(zeros(m,n));
xbar = zeros(m,n);

tic;
for it = 1:iter_max    
    xi = prox_fstar(sigma, xi +sigma*A(xbar));
    xold = x;
    x = prox_g(x-tau*Astar(xi));
    xbar = x+theta*(x-xold);
end
computing_time = toc;
function [x, computing_time, l1values] = DouglasRachfordMin(x0, prox_model, prox_fidelity, options)
% x0 is the initial point
% prox_model is the proximal to the data model (say l1)
% prox_fidelity is the projection operator on the data fidelity (say subsampled Fourier transform of x is y)
% iter_max: the maximum number of iterations 
% lambda: the relaxation parameter
% gamma: the 'order' of the method
%
% This function minimizes model(x) + fidelity(x) with the following
% updates: 
% 1. u(k+1) = prox_fidelity_gamma(x(k))
% 2. v(k+1) = prox_model_gamma(2u(k+1)-x(k)
% 3. x(k+1) = x(k) + lambda(v(k+1)-u(k+1))

%% Check the inputs 

if isempty(x0)
    disp("DouglasRachfordMin:x0 is empty ... program about to crash!")
end


iter_max = 200;
lambda = 1.5;
gamma = 0.05;
verbose = 20;


if ~(nargin < 4 || isempty(options))
    if isfield(options, 'iter_max')
        iter_max = options.iter_max;
    end
    if isfield(options, 'lambda')
        lambda = options.lambda;
    end
    if isfield(options, 'gamma')
        gamma = options.gamma;
    end
    if isfield(options, 'verbose')
        verbose = options.verbose;
    end
end

% Random things to set
l1values = zeros(iter_max);
x = x0; 

%% Main iterations 
tic;
for cur_iter=1:iter_max % Need to add a stopping criterion based on dual gap or something like this. 
    % Step 1: check the data fidelity: 
    u = prox_fidelity(x,gamma);
    % Step 2: enforce the model consistency: 
    v = prox_model(2*u-x,gamma);
    % Step 3: "combine by regularization":
    x = x + lambda *(v-u);
    
    l1values(cur_iter) = sum(abs(x(:)));
%     imagesc([u,v,x])
%     colormap(gray(256))
%     pause(3)
    
    if mod(cur_iter, verbose) == 0
        disp(['Current iteration: ', num2str(cur_iter), '. Norm of the current estimate is ', num2str(l1values(cur_iter))]);
        
    end

end
computing_time = toc;
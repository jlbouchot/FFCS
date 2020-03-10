function [x_hat] = CS_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = zeros(size(measurementStar(y_rhs)));
end
disp('Recovery by compressed sensing method')
% In traditional compressed sensing, there are no operators inside the
% first norm. Hence, A and Astar are simply the identity. 
A = @(f) f;
Astar = @(f) f; 


% Proximal operator for the l1 norm
prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
% Proximal operator for the l2 data fidelity term
prox_G = @(z) z + real(measurementStar(y_rhs- f_mask.*measurement(z)));


% Call the primal dual algorithm
disp('          ... primal dual algorithm running ... ')

x_hat = PrimalDualMin(x0, prox_Fstar, prox_G, A, Astar, PDoptions);


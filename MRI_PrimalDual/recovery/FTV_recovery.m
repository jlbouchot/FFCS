function [x_hat] = FTV_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = zeros(size(measurementStar(y_rhs{1})));
end
disp('Recovery by compressed sensing method')
% In traditional compressed sensing, there are no operators inside the
% first norm. Hence, A and Astar are simply the identity. 


[m,n] = size(x0);

A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));

abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
% Proximal operator for the l1 norm
prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

nb_channels = length(y_rhs);
for one_channel=1:nb_channels
    
    % Proximal operator for the l2 data fidelity term
    prox_G = @(z) z + real(measurementStar(y_rhs{one_channel}- f_mask.*measurement(z)));
    
    
    % Call the primal dual algorithm
    disp('          ... primal dual algorithm running ... ')
    
    x_hat{one_channel} = PrimalDualMin(x0, prox_Fstar, prox_G, A, Astar, PDoptions);
end

function f_mask = sample_Fourier_VDS(N, deterministic_ratio)
% Samples a deterministic core region. 
% This is related to some work by Chauffert et all  on traveling salesman 
% optimization and variable density sampling in MRI reconstruction. 

if nargin < 2 || isempty(deterministic_ratio)
    deterministic_ratio = 1/3; % in one direction. does not correpond to ratio in the 2D space.
end

X = ones(N,1)*linspace(-1,1,N);
Y = linspace(-1,1,N)'*ones(1,N);
pts_distances = X.^2+Y.^2;

f_mask = (pts_distances < deterministic_ratio);
return
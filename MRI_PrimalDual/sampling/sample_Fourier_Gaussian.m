function f_mask = sample_Fourier_Gaussian(gaussian_ratio, N)

X = ones(N,1)*linspace(-1,1,N);
Y = linspace(-1,1,N)'*ones(1,N);
rejection_threshold = exp(-gaussian_ratio*X.^2 - gaussian_ratio*Y.^2);
rejection_threshold = rejection_threshold / max(max(rejection_threshold));

f_mask = (rand(N,N) < rejection_threshold); % This generate a mask where most of the non zeros are at the (0,0) frequencies
return
function f_mask = sample_Fourier_uniform(rejection_threshold, N)

f_mask = (rand(N,N) <= rejection_threshold);

return
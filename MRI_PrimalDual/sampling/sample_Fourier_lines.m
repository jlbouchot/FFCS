function f_mask = sample_Fourier_lines(nb_line_per_dim, size_per_dim)

dim1_samples = min(max(unique(sort(floor(sqrt(1.5*size_per_dim)*randn(nb_line_per_dim,1))))+size_per_dim/2, 0), size_per_dim);
dim2_samples = min(max(unique(sort(floor(sqrt(1.5*size_per_dim)*randn(nb_line_per_dim,1))))+size_per_dim/2, 0), size_per_dim);

f_mask = zeros(size_per_dim);
f_mask(dim1_samples,:) = 1;
f_mask(:,dim2_samples) = 1; 

return
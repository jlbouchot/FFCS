function f_mask = sample_Fourier_oneD(nb_lines, size_per_dim, go_vertical)

sigma = sqrt(size_per_dim);

% Sample without replacement
values = 1:size_per_dim; 
Rweights = zeros(size_per_dim,1);
% for one_idx=1:size_per_dim
%     k = one_idx-1;
%     weights(one_idx) = (target/p)^k*(1-target/p)^(size_per_dim-1-k)*nchoosek(size_per_dim-1,k);
% end
weights = exp(-(values-(size_per_dim+1)/2).^2/sigma^2);


% Use "ceil" below to make sure we have an integer!
selected_idx = datasample(values, ceil(nb_lines), 'Replace', false, 'Weights', weights);

% Create the mask
f_mask = zeros(size_per_dim);

if go_vertical
    f_mask(:,selected_idx) = 1;
else
    f_mask(selected_idx,:) = 1;
end

return
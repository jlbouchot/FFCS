function rec_proj = proj_on_details(x_i, options)

nb_lvl = options.wLevel;
cur_lvl = options.detailLevel;
wname = options.wname;

% Have to add the handling of the various parameters

[wave_coefs,bookkeeping] = wavedec2(x_i,nb_lvl,wname);

tot_scales = size(bookkeeping,1)-2;
dim_per_scale = prod(bookkeeping,2);
cur_details = zeros(size(wave_coefs));
cur_dim = dim_per_scale(end-cur_lvl,:);

if tot_scales == cur_lvl
    start_idx = 1;
    cur_details(start_idx:4*cur_dim) = wave_coefs(start_idx:4*cur_dim);
else
    start_idx = dim_per_scale(1) + 3*sum(dim_per_scale(2:end-1-cur_lvl));
    cur_details([start_idx+1:start_idx+3*cur_dim]) = wave_coefs([start_idx+1:start_idx+3*cur_dim]);
    % cur_details([1:dim_per_scale(1), start_idx+1:start_idx+3*cur_dim]) = wave_coefs([1:dim_per_scale(1), start_idx+1:start_idx+3*cur_dim]);
end

rec_proj = waverec2(cur_details,bookkeeping,wname);

end
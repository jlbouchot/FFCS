function filtered_x = canonical_projection(x,options)

proj_mask = options.mask;
filtered_x = x.*proj_mask;

return
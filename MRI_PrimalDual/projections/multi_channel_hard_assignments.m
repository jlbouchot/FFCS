function light_masks = multi_channel_hard_assignments(N)
% For now, we only deal with 6 channels: 
% Center quarter and ring 3/4, 4 corner quarters

a_quarter = floor(N/4);
a_half = floor(N/2);

for one_light=1:6
    light_masks{one_light} = zeros(N,N);
end

light_masks{1}(a_quarter+1:a_half+a_quarter,a_quarter+1:a_half+a_quarter) = 1; 
light_masks{2} = 1 - light_masks{1};
light_masks{3}(1:a_half,1:a_half) = 1; 
light_masks{4}(1:a_half,a_half+1:end) = 1; 
light_masks{5}(1+a_half:end,a_half+1:end) = 1; 
light_masks{6}(1+a_half:end,1:a_half) = 1; 

return
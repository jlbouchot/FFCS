function light_masks = multi_channel_sines(N, nb_sines)

if nargin < 2 || isempty(nb_sines)
    nb_sines = 5;
end

start_x = ceil(N/(nb_sines)*(0:nb_sines))+1;
max_idx = max(max(start_x(3:end) - start_x(1:end-2)),start_x(2) + start_x(end)-start_x(end-1));
x_values = 1:max_idx; 

cosine_values = 1- cos(nb_sines*pi*x_values/N);

for one_sine = 1:nb_sines-1
    light_masks{one_sine} = zeros(N);
    cur_x = start_x(one_sine):start_x(one_sine+2)-1;
    light_masks{one_sine}(:,cur_x) = ones(N,1)*cosine_values(1:start_x(one_sine+2)-start_x(one_sine));
end

light_masks{nb_sines} = zeros(N,N);
nb_idx_left = N- start_x(end-1)+1;
light_masks{nb_sines}(:,start_x(end-1):end) = ones(N,1)*cosine_values(1:nb_idx_left);
light_masks{nb_sines}(:,1:start_x(2)-1) = ones(N,1)*cosine_values(nb_idx_left+1:nb_idx_left+start_x(2)-1);


return
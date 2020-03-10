function light_masks = multi_channel_lightning(N, nb_lights, centers, variance)

if nargin < 2 || isempty(nb_lights)
    nb_lights = 4;
    centers = [floor(N/4), floor(N/4); floor(N/4), floor(3*N/4); floor(3*N/4), floor(N/4); floor(3*N/4), floor(3*N/4)];
end

if nargin < 3 || isempty(centers)
   centers = floor(N*rand(nb_lights,2)); 
end

if nargin < 4 || isempty(variance)
    variance = N*N;
end

for one_light=1:nb_lights
    X = ones(N,1)*((1:N)-centers(one_light,1));
    Y = ((1:N)- centers(one_light,2))'*ones(1,N);
    light_masks{one_light} = exp(-(X.^2+Y.^2)/variance);
end

return
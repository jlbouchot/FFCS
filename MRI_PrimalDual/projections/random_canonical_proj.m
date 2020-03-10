function [masks, S] = random_canonical_proj(nb_proj, N, rank)

S = zeros(N); % Will contain the frame operator

for one_mask = 1:nb_proj
    masks{one_mask} = (rand(N) <= rank/N/N);
    S = S + masks{one_mask};
end

return
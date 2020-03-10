function masks = ring_partitions(nb_rings, N, type)

% type will be either 'uniform' or 'growing'
% The uniform pretty much keeps the same radius for all rings up to the
% last one which contains everything that is left.
% The growing one will contain rings of width 1r_0, 2r_0, ... nb-1 r_0
% until the last one which will contain the corners. 

if mod(N,2) == 0
    X = ones(N,1)*(-N/2+1:N/2);
    Y = (-N/2+1:N/2)'*ones(1,N);
else
    X = ones(N,1)*(-N/2:N/2);
    Y = (-N/2:N/2)'*ones(1,N);
end

Z = sqrt(X.^2 + Y.^2);

masks = ones(N);

if strcmp(type,'uniform')
    splits = 0:N/2/(nb_rings-1):N/2;
%     for one_ring=nb_rings-1:-1:1
%         sum(sum(Z < splits(one_ring+1)))
%         Z(Z < splits(one_ring+1)) = one_ring+1;
%     end
    
    for one_ring=1:nb_rings
%         sum(sum(Z >= splits(one_ring)))
        masks(Z >= splits(one_ring)) = one_ring;
    end
else
    
end


return
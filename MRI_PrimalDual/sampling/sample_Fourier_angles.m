function f_mask = sample_Fourier_angles(nb_lines, N)

rdm_angles = 2*pi*rand(nb_lines,1);

f_mask = zeros(N);
if mod(N,2) == 0
    t = -N/2+1:N/2;
    t_shift = N/2;
else
    t = -N/2:N/2;
    t_shift = N/2+1; % This will allow to reset the right coordinate. 
end

for one_angle=1:nb_lines
    alpha = rdm_angles(one_angle);
    if alpha < pi/4 || (alpha >= 3*pi/4 && alpha < 5*pi / 4) || alpha >= 7*pi/4
        x = t; 
        y = round(t*tan(alpha));
    else % Technically, we should distinguishe the case alpha = +/- pi/2 -- which will happen with probability 0. 
        y = t; 
        x = round(t/tan(alpha));
    end
    for one_point = 1:length(t)
        f_mask(y(one_point)+t_shift,x(one_point)+t_shift) = 1;
    end
end

return
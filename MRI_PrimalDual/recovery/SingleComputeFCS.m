function x_hat = SingleComputeFCS(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions, masks)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = ones(size(measurementStar(y_channel{1})));
end
disp('Single sensor Fused L1 analysis recovery')

if nargin < 7 || isempty(WAVEoptions)
    wave_name = WAVEoptions.wname;
    nb_level = WAVEoptions.nblvl;
    
    % Get the structure first, this will not change through time
    [~,s]=wavedec2(x0,nb_level,wave_name);
    
    % Analysis operator
    A = @(f) wavedec2(f,nb_level,wave_name);
    % Synthesis operator
    Astar = @(f) waverec2(f,s,wave_name);
else
    A = @(f) f;
    Astar = @(f) f;
end

nb_channels = length(y_channel);
y_rhs = [];
S = zeros(size(masks{1}));
for one_channel=1:nb_channels
    y_rhs = [y_rhs; y_channel{one_channel}];
    S = S+masks{one_channel};
end


% Proximal operator for the l1 norm
prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
% Proximal operator for the l2 data fidelity term
% prox_G = @(z) z + real(measurementStar(y_rhs- multichannelStacking(Astar, f_mask, masks, z)));
prox_G = @(z) z + real(multichannelUnstacking(y_rhs- multichannelStacking(Astar, f_mask, masks, z), ... 
    nb_channels, A, S));
% multichannelUnstacking(input_y, nb_channels, A, S)

x_hat = Astar(PrimalDualMin(x0, prox_Fstar, prox_G, @(f) f, @(f) f, PDoptions));

end

function output = multichannelStacking(WS_op, f_mask, masks, x)
nb_masks = length(masks); 

output = [];
synthesised_x = WS_op(x);

for one_channel=1:nb_masks
    output = [output;f_mask.*c_fft_2d( masks{one_channel}.*synthesised_x)];
end

end

function output = multichannelUnstacking(input_y, nb_channels, WA_op, S)

output_size = [size(input_y,1)/nb_channels, size(input_y,2)];


values = zeros(output_size(1),output_size(2), nb_channels);
for one_channel = 1:nb_channels
    values(:,:,one_channel) = c_ifft_2d(input_y((one_channel-1)*output_size(1)+1:one_channel*output_size(1),:));
end

output = WA_op(sum(values,3)./S);


end


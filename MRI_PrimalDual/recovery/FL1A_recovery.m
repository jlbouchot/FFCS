function x_hat = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = ones(size(measurementStar(y_channel{1})));
end
disp('Fused L1 analysis recovery')

if nargin < 7 || isempty(WAVEoptions)
    A = @(f) f;
    Astar = @(f) f;
else
    wave_name = WAVEoptions.wname;
    nb_level = WAVEoptions.nblvl;
    
    % Get the structure first, this will not change through time
    [~,s]=wavedec2(x0,nb_level,wave_name);
    
    % Analysis operator
    A = @(f) wavedec2(f,nb_level,wave_name);
    % Synthesis operator
    Astar = @(f) waverec2(f,s,wave_name);
end

prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

nb_channels = length(y_channel);

for one_channel=1:nb_channels
    prox_G = @(z) z + real(measurementStar(y_channel{one_channel}- f_mask.*measurement(z)));
%     prox_G = @(z) z + real(c_ifft_2d(y_channel{one_channel}- f_mask.*c_fft_2d(z)));
    
    % Call the primal dual algorithm
    disp('          ... primal dual algorithm running ... ')
    x_hat{one_channel} = PrimalDualMin(x0, prox_Fstar, prox_G, A, Astar, PDoptions);
end
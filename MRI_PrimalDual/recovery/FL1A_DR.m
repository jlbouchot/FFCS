function x_local = FL1A_DR(y_channel, f_mask, x0, measurement, measurementStar, WS_op, WA_op, DRoptions)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = ones(size(measurementStar(y_channel{1})));
end
disp("Computing the Fused l1 analysis recovery with Douglas-Rachford.. ")

nb_channels = length(y_channel);

for one_channel=1:nb_channels
    prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    prox_l1_analysis = @(z,sigma) WS_op(prox_l1_straight(WA_op(z),sigma));
    
    A=@(f) f_mask.*(measurement(f));
    At = @(f,gamma) f +real(measurementStar(y_channel{one_channel}-A(f)));
    %         At = @(f,gamma) f +c_ifft_2d(y_channel{one_channel}-A(f));
    
    
    % Call the Douglas-Rachford algorithm
    % disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
    x_local{one_channel} = DouglasRachfordMin(x0,prox_l1_analysis,At,DRoptions);
end



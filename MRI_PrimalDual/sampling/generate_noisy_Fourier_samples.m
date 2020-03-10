function y_out = generate_noisy_Fourier_samples(img, noise_lvl, noiseName, f_mask, data_masks)

if nargin < 5 || isempty(data_masks)
    data_masks{1} = ones(size(img));
end

nb_channels = length(data_masks);
[m,n] = size(img);
for one_channel=1:nb_channels
    switch noiseName
        case 'Normal'
            img_f = c_fft_2d(data_masks{one_channel}.*img) + noise_lvl*randn(m,n);
        case 'Speckle'
            img_f = imnoise(c_fft_2d(data_masks{one_channel}.*img), 'speckle', noise_lvl);
        case 'SnP'
            img_f = imnoise(c_fft_2d(data_masks{one_channel}.*img), 'salt & pepper', noise_lvl);
    end
    y_out{one_channel} = f_mask.*img_f;
    
end
if nb_channels == 1 
    y_out = y_out{1};
end
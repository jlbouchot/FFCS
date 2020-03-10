function band_filtered_x = proj_Fourier_band(x, options)

f_mask = options.mask;
band_filtered_x = c_ifft_2d(f_mask.*c_fft_2d(x));

return
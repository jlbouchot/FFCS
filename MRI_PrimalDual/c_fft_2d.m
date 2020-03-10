function f_of_xi = c_fft_2d(X) % Stands for "centered FFT in 2D" 
f_of_xi = 1/sqrt(numel(X))*fftshift(fft2(ifftshift(X)));
end

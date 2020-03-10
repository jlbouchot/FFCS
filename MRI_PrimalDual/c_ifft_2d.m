function f_of_x = c_ifft_2d(X_hat)
f_of_x = sqrt(numel(X_hat))*ifftshift(ifft2(fftshift(X_hat)));
end
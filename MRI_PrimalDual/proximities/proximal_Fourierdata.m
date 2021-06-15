function [prox_FourierData, prox_FourierDataStar] = proximal_Fourierdata(z, tau)

prox_FourierData = (real(c_ifft_2d(y_rhs)) + z) - ... 
    tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_rhs) + z)));
prox_FourierDataStar = z-tau*prox_FourierData(1/tau, z/tau);

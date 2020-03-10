% This script calls the recovery from the variable density paper. 
% 
% It includes the following sampling patterns 
% * VDS with Gaussian samples 
% * VDS with a spiral sample 
% 
% Recovery using 
% * DR from Chauffert 
% * DR for fused L1 analysis 
% * PD for fused TV minimization 
% 
% The tests are made on 2d images such as 
% * The brain used in Chauffert et al 
% * The brain from Michi Lustig 
% 
% As a result, we have 
% * All the recovery saved 
% * All the sampling patterns saved 
% * All the pointwise errors saved 
% * All the quality metrics (PSNR, SSIM, L2) measured 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Bits and pieces
clear all
close all
clc

addpath(genpath('recovery/')) % Contains the self implemented recovery
addpath(genpath('data/')) % Contains the brain images
addpath(genpath('VDS_Chauffert/')) % Contains the DR L1 min and sampling
addpath(genpath('projections/')) % Used for the fused methods
addpath('./utils/') 

% Toolboxes used by Chauffert 
addpath(genpath('toolbox_general/'))
addpath(genpath('toolbox_signal/'))

% Let's go 
disp('Using variable density for MR Image denoising')


%% Parameters for both sampling scheme
load reference      % load I, a reference image
N_I=size(I);
I = (I - min(I(:)))/(max(I(:)) - min(I(:)));

load brain          % Loads im, another reference image
N_im = size(im);
im = abs(im);
im = (im - min(im(:)))/(max(im(:)) - min(im(:)));


% Specify wavelet options
options.h= compute_wavelet_filter('Symmlet',10);
myWT=@(x) perform_wavortho_transf(x,4,1,options);
myIWT=@(x) perform_wavortho_transf(x,4,-1,options);

%% Choose sampling scheme parameters:

% Target distribution
opts1.distrib=compute_optimal_distrib(N_I,options,[3 3]);     % Computes optimal distribution
opts2.distrib=compute_optimal_distrib(N_im,options,[3 3]);     

% A deterministic part:
opts1.deterministic_part=set_LF(opts1.distrib,3);
opts2.deterministic_part=set_LF(opts2.distrib,3);


% Choose the type of sampling scheme, could be:
% 'TSP', 'Independent', 'Radial', 'Radial_Random',
% 'Markov' or 'Spiral'.
opts1.type='Radial Random';
opts2.type = 'Radial Random';

opts1spiral = opts1;
opts1spiral.type = 'Spiral';

opts2spiral = opts2; 
opts2spiral.type = 'Spiral';
% Sampling ratio 
R=10;


%% Generate sampling schemes
samplingScheme1I=generate2D_CS_scheme(N_I,1/R,opts1);
f_mask_I1 = samplingScheme1I.mask;
samplingScheme2I=generate2D_CS_scheme(N_I,1/R,opts1spiral);
f_mask_I2 = samplingScheme2I.mask;

samplingScheme1im=generate2D_CS_scheme(N_im,1/R,opts2);
f_mask_im1 = samplingScheme1im.mask;
samplingScheme2im=generate2D_CS_scheme(N_im,1/R,opts2spiral);
f_mask_im2 = samplingScheme2im.mask;

%% Generate the prefiltering operations
nb_channels = 8;
masks_I = multi_channel_sines(N_I(1), nb_channels);
S_I = zeros(N_I);
for one_channel=1:nb_channels
    S_I = S_I + masks_I{one_channel};
end

masks_im = multi_channel_sines(N_im(1), nb_channels);
S_im = zeros(N_im);
for one_channel=1:nb_channels
    S_im = S_im + masks_im{one_channel};
end

%% Generate measurements
noise_lvl = 0.05;
noiseName = 'Normal';
% First for the multi-sensor case
y_channel_I1 = generate_noisy_Fourier_samples(I, noise_lvl, noiseName, f_mask_I1, masks_I);
y_channel_im1 = generate_noisy_Fourier_samples(im, noise_lvl, noiseName, f_mask_im1, masks_im);
% Then for the single sensor case
y_rhs_I1 = generate_noisy_Fourier_samples(I, noise_lvl, noiseName, f_mask_I1);
y_rhs_im1 = generate_noisy_Fourier_samples(im, noise_lvl, noiseName, f_mask_im1);


% First for the multi-sensor case
y_channel_I2 = generate_noisy_Fourier_samples(I, noise_lvl, noiseName, f_mask_I2, masks_I);
y_channel_im2 = generate_noisy_Fourier_samples(im, noise_lvl, noiseName, f_mask_im2, masks_im);
% Then for the single sensor case
y_rhs_I2 = generate_noisy_Fourier_samples(I, noise_lvl, noiseName, f_mask_I2);
y_rhs_im2 = generate_noisy_Fourier_samples(im, noise_lvl, noiseName, f_mask_im2);


%% Recovery via Chauffert 
disp('-----------------------')
disp('Reconstructing MR image')
disp('-----------------------')

% The first one
A=@(x) f_mask_I1.*(c_fft_2d(myIWT(x)));
At=@(x) myWT(c_ifft_2d(x));

I_chauffert1 = abs(myIWT(Solve_l1_problemDR(y_rhs_I1,A,At,200)));

A=@(x) f_mask_I2.*(c_fft_2d(myIWT(x)));
At=@(x) myWT(c_ifft_2d(x));
I_chauffert2 = abs(myIWT(Solve_l1_problemDR(y_rhs_I2,A,At,200)));

% The second one 
A=@(x) f_mask_im1.*(c_fft_2d(myIWT(x)));
At=@(x) myWT(c_ifft_2d(x));

im_chauffert1 = abs(myIWT(Solve_l1_problemDR(y_rhs_im1,A,At,200)));

A=@(x) f_mask_im2.*(c_fft_2d(myIWT(x)));
At=@(x) myWT(c_ifft_2d(x));

im_chauffert2 = abs(myIWT(Solve_l1_problemDR(y_rhs_im2,A,At,200)));

%% Recovery with fused l1analysis

options.wname = 'Symmlet';
options.h= compute_wavelet_filter(options.wname,10);
WA_op=@(x) perform_wavortho_transf(x,4,1,options);
WS_op=@(x) perform_wavortho_transf(x,4,-1,options);
DRoptions.gamma = 0.05; % Used for the soft-thresholding part of the update
DRoptions.lambda =1.5; % Used for the relaxation of the update
DRoptions.iter_max = 180;
DRoptions.verbose = 30;

% Radial lines plus low-frequency sampling
x_local = FL1A_DR(y_channel_I1, f_mask_I1, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), WS_op, WA_op, DRoptions);
x_fl1a_I1 = zeros(N_I);
for one_channel=1:nb_channels
    x_fl1a_I1 = x_fl1a_I1 + x_local{one_channel};
end
x_fl1a_I1 = abs(x_fl1a_I1./S_I);


x_local = FL1A_DR(y_channel_im1, f_mask_im1, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), WS_op, WA_op, DRoptions);
x_fl1a_im1 = zeros(N_im);
for one_channel=1:nb_channels
    x_fl1a_im1 = x_fl1a_im1 + x_local{one_channel};
end
x_fl1a_im1 = abs(x_fl1a_im1./S_im);

% Spiral with low frequency sampling
x_local = FL1A_DR(y_channel_I2, f_mask_I2, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), WS_op, WA_op, DRoptions);
x_fl1a_I2 = zeros(N_I);
for one_channel=1:nb_channels
    x_fl1a_I2 = x_fl1a_I2 + x_local{one_channel};
end
x_fl1a_I2 = abs(x_fl1a_I2./S_I);


x_local = FL1A_DR(y_channel_im2, f_mask_im2, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), WS_op, WA_op, DRoptions);
x_fl1a_im2 = zeros(N_im);
for one_channel=1:nb_channels
    x_fl1a_im2 = x_fl1a_im2 + x_local{one_channel};
end
x_fl1a_im2 = abs(x_fl1a_im2./S_im);


% figure
% colormap gray
% subplot(1,3,1), imagesc(f_mask_I1), axis equal,axis off
% title(['mask: ' opts1.type])
% subplot(1,3,2), imagesc(x_fl1a_I1), axis equal, axis off
% title('reconstructed image')
% subplot(1,3,3), imagesc(abs(I-x_fl1a_I1)), axis equal, axis off
% title(['difference with original: (psnr= ' num2str(psnr(x_fl1a_I1,I)) ')'])
% 
% figure
% colormap gray
% subplot(1,3,1), imagesc(f_mask_I1), axis equal,axis off
% title(['mask: ' opts1.type])
% subplot(1,3,2), imagesc(I_chauffert1), axis equal, axis off
% title('reconstructed image')
% subplot(1,3,3), imagesc(abs(I-I_chauffert1)), axis equal, axis off
% title(['difference with original: (psnr= ' num2str(psnr(I_chauffert1,I)) ')'])
%% Recovery with fused TV 
PDoptions.iter_max = 200;
PDoptions.theta = 1;
PDoptions.sigma = 2;
PDoptions.tau = .95/16/PDoptions.sigma;

% Radial lines and low frequency samplimg
x_local = FTV_recovery(y_channel_I1,f_mask_I1, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), PDoptions);
x_TV_I1 = zeros(N_I);
for one_channel=1:nb_channels
    x_TV_I1 = x_TV_I1 + x_local{one_channel};
end
x_TV_I1 = abs(x_TV_I1./S_I);

x_local = FTV_recovery(y_channel_im1,f_mask_im1, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), PDoptions);
x_TV_im1 = zeros(N_im);
for one_channel=1:nb_channels
    x_TV_im1 = x_TV_im1 + x_local{one_channel};
end
x_TV_im1 = abs(x_TV_im1./S_im);



% Spiral and low frequency samplimg
x_local = FTV_recovery(y_channel_I2,f_mask_I2, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), PDoptions);
x_TV_I2 = zeros(N_I);
for one_channel=1:nb_channels
    x_TV_I2 = x_TV_I2 + x_local{one_channel};
end
x_TV_I2 = abs(x_TV_I2./S_I);

x_local = FTV_recovery(y_channel_im2,f_mask_im2, [], @(f) c_fft_2d(f), @(f) c_ifft_2d(f), PDoptions);
x_TV_im2 = zeros(N_im);
for one_channel=1:nb_channels
    x_TV_im2 = x_TV_im2 + x_local{one_channel};
end
x_TV_im2 = abs(x_TV_im2./S_im);




% 
% figure
% colormap gray
% subplot(1,3,1), imagesc(f_mask_I1), axis equal,axis off
% title(['mask: ' opts1.type])
% subplot(1,3,2), imagesc(x_TV_I1), axis equal, axis off
% title('reconstructed image')
% subplot(1,3,3), imagesc(abs(I-x_TV_I1)), axis equal, axis off
% title(['difference with original: (psnr= ' num2str(psnr(x_TV_I1,I)) ')'])

%% Displays 
% Chauffert's method
figure('name', 'Chauffert - Radial')
colormap gray
subplot(2,3,1), imagesc(f_mask_I1), axis equal,axis off
title(['mask: ' opts1.type])
subplot(2,3,2), imagesc(I_chauffert1), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-I_chauffert1)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(I_chauffert1,I)) ', SSIM = ', num2str(ssim(I_chauffert1,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im1), axis equal,axis off
title(['mask: ' opts2.type])
subplot(2,3,5), imagesc(im_chauffert1), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-im_chauffert1)), axis equal, axis off
% title(['difference with original: (psnr= ' num2str(psnr(im_chauffert1,im)) ')'])
title(['Error: (psnr= ' num2str(psnr(im_chauffert1,im)) ', SSIM = ', num2str(ssim(im_chauffert1,im)), ')'])

figure('name', 'Chauffert - Spiral')
colormap gray
subplot(2,3,1), imagesc(f_mask_I2), axis equal,axis off
title(['mask: ' opts1spiral.type])
subplot(2,3,2), imagesc(I_chauffert2), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-I_chauffert2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(I_chauffert2,I)) ', SSIM = ', num2str(ssim(I_chauffert2,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im2), axis equal,axis off
title(['mask: ' opts2spiral.type])
subplot(2,3,5), imagesc(im_chauffert2), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-im_chauffert2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(im_chauffert2,im)) ', SSIM = ', num2str(ssim(im_chauffert2,im)), ')'])

% Fused l1 analysis
figure('name', 'Fused l1 analysis - Radial')
colormap gray
subplot(2,3,1), imagesc(f_mask_I1), axis equal,axis off
title(['mask: ' opts1.type])
subplot(2,3,2), imagesc(x_fl1a_I1), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-x_fl1a_I1)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_fl1a_I1,I)) ', SSIM = ', num2str(ssim(x_fl1a_I1,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im1), axis equal,axis off
title(['mask: ' opts2.type])
subplot(2,3,5), imagesc(x_fl1a_im1), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-x_fl1a_im1)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_fl1a_im1,im)) ', SSIM = ', num2str(ssim(x_fl1a_im1,im)), ')'])


figure('name', 'Fused l1 analysis - Spiral')
colormap gray
subplot(2,3,1), imagesc(f_mask_I2), axis equal,axis off
title(['mask: ' opts1.type])
subplot(2,3,2), imagesc(x_fl1a_I2), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-x_fl1a_I2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_fl1a_I2,I)) ', SSIM = ', num2str(ssim(x_fl1a_I2,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im2), axis equal,axis off
title(['mask: ' opts2.type])
subplot(2,3,5), imagesc(x_fl1a_im2), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-x_fl1a_im2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_fl1a_im2,im)) ', SSIM = ', num2str(ssim(x_fl1a_im2,im)), ')'])


% Fused TV analysis
figure('name', 'Fused TV minimization - Radial')
colormap gray
subplot(2,3,1), imagesc(f_mask_I1), axis equal,axis off
title(['mask: ' opts1.type])
subplot(2,3,2), imagesc(x_TV_I1), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-x_TV_I1)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_TV_I1,I)) ', SSIM = ', num2str(ssim(x_TV_I1,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im1), axis equal,axis off
title(['mask: ' opts2.type])
subplot(2,3,5), imagesc(x_TV_im1), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-x_TV_im1)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_TV_im1,im)) ', SSIM = ', num2str(ssim(x_TV_im1,im)), ')'])


figure('name', 'Fused TV minimization - Spiral')
colormap gray
subplot(2,3,1), imagesc(f_mask_I2), axis equal,axis off
title(['mask: ' opts1.type])
subplot(2,3,2), imagesc(x_TV_I2), axis equal, axis off
title('reconstructed image')
subplot(2,3,3), imagesc(abs(I-x_TV_I2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_TV_I2,I)) ', SSIM = ', num2str(ssim(x_TV_I2,I)), ')'])
subplot(2,3,4), imagesc(f_mask_im2), axis equal,axis off
title(['mask: ' opts2.type])
subplot(2,3,5), imagesc(x_TV_im2), axis equal, axis off
title('reconstructed image')
subplot(2,3,6), imagesc(abs(im-x_TV_im2)), axis equal, axis off
title(['Error: (psnr= ' num2str(psnr(x_TV_im2,im)) ', SSIM = ', num2str(ssim(x_TV_im2,im)), ')'])

%% Save everything

extensions = {'jpg', 'png', 'tif'};


vdsOutput_path = fullfile('.', 'results', 'VDSComparison');

if ~isfolder(vdsOutput_path)
    mkdir(vdsOutput_path);
end
 

for one_extension=1:length(extensions)
    % Start with the cameraman experiment -> camMan_dir_path
    % Save the sampling patterns
    curRes = f_mask_I1; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_CScan_Radial.', extensions{one_extension}]))
    curRes = f_mask_I2; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_CScan_Spiral.', extensions{one_extension}]))
    curRes = f_mask_im1; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_RealBrain_Radial.', extensions{one_extension}]))
    curRes = f_mask_im2; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_RealBrain_Spiral.', extensions{one_extension}]))
    
    % Save the original images 
    curRes = I; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_Original.', extensions{one_extension}]))
    curRes = im; 
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_Original.', extensions{one_extension}]))
    
    % Save the recoveries 
    % and the pointwise errors
    % Original work from Chauffert
    curRes = abs(I_chauffert1); 
    curErr = abs(I_chauffert1 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_VDS_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_VDS_Radial_error.', extensions{one_extension}]))
    curRes = abs(I_chauffert2); 
    curErr = abs(I_chauffert2 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_VDS_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_VDS_Spiral_error.', extensions{one_extension}]))
    
    curRes = abs(im_chauffert1); 
    curErr = abs(im_chauffert1 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_VDS_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_VDS_Radial_error.', extensions{one_extension}]))
    curRes = abs(im_chauffert2); 
    curErr = abs(im_chauffert2 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_VDS_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_VDS_Spiral_error.', extensions{one_extension}]))
    
    % Fused L1 analysis
    curRes = abs(x_fl1a_I1); 
    curErr = abs(x_fl1a_I1 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_FL1_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_FL1_Radial_error.', extensions{one_extension}]))
    curRes = abs(x_fl1a_I2); 
    curErr = abs(x_fl1a_I2 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_FL1_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_FL1_Spiral_error.', extensions{one_extension}]))
    
    curRes = abs(x_fl1a_im1); 
    curErr = abs(x_fl1a_im1 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_FL1_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_FL1_Radial_error.', extensions{one_extension}]))
    curRes = abs(x_fl1a_im2); 
    curErr = abs(x_fl1a_im2 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_FL1_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_FL1_Spiral_error.', extensions{one_extension}]))
    
    % Fused TV analysis
    curRes = abs(x_TV_I1); 
    curErr = abs(x_TV_I1 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_FTV_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_FTV_Radial_error.', extensions{one_extension}]))
    curRes = abs(x_TV_I2); 
    curErr = abs(x_TV_I2 - I);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_FTV_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['CScan_FTV_Spiral_error.', extensions{one_extension}]))
    
    curRes = abs(x_TV_im1); 
    curErr = abs(x_TV_im1 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_FTV_Radial_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_FTV_Radial_error.', extensions{one_extension}]))
    curRes = abs(x_TV_im2); 
    curErr = abs(x_TV_im2 - im);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_FTV_Spiral_Recovery.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['RealBrain_FTV_Spiral_error.', extensions{one_extension}]))
   
end


save(fullfile(vdsOutput_path, 'allResults.mat'))
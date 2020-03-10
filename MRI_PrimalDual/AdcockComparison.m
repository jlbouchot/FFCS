% Tests used in the paper:
%
%
addpath('./sampling/') % Contains the various sampling schemes
addpath('./projections/') % Contains various filtering operations
addpath('./data/') % Some examples of data
addpath('./recovery/') % Some recovery procedures
addpath('./utils/') % Some recovery procedures

%% Some cleaning
clear all
close all
clc

%% Parameters for all to share

% Parameters for the Primal dual algorithm
% These values do not change throughout the experiments.
PDoptions.iter_max = 200;
PDoptions.theta = 1;
PDoptions.sigma = 2;
PDoptions.tau = .95/16/PDoptions.sigma;

tau = PDoptions.tau;

% Wavelet basis
WAVEoptions.wname = 'db4';
WAVEoptions.nblvl = 3;

%% Which test to run?
do_cameraman = true; % Run the tests for the natural scene

%% First test: Natural scene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This experiment handles the cameraman image and checks its recovery from
% multiple sensors / subsampled images.
%
% The recovery is done using a primal dual approach and proximal operators.
% The sampling is done by subsampling the Fourier transform of the image
% according to a Gaussian distribution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cameraman image
if do_cameraman
    % Load the image
    img_cam = double(imread('cameraman.tif'));
    [m,n] = size(img_cam);
    img_cam = img_cam/max(img_cam(:));
    
    % Generate sampling pattern
    how_spread = m^(1/2); % The bigger this parameter, the fewer sampling points will be taken
    f_mask = sample_Fourier_Gaussian(how_spread, m);
    
    % Generate sensor pre-filtering
    nb_channels = 8;
    masks = multi_channel_sines(m, nb_channels);
    S = zeros(size(img_cam));
    for one_channel=1:nb_channels
        S = S + masks{one_channel};
    end
    
    % Generate measurements
    noise_lvl = 0.01;
    noiseName = 'Speckle';
    % First for the multi-sensor case
    y_channel = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask, masks);
    % Then for the single sensor case
    y_rhs = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask);
    
    % Now we go ahead with the actual recovery
    x0 = [];
    measurement = @(f) c_fft_2d(f);
    measurementStar = @(f) c_ifft_2d(f);
    x_adcock = SingleComputeFCS(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions, masks);
    x_cs = CS_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions);
    x_l1a = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    x_fl1a_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    % Compute the fused l1 analysis solution:
    x_fl1a = zeros(size(img_cam));
    for one_channel=1:nb_channels
        x_fl1a = x_fl1a + x_fl1a_all{one_channel};
    end
    x_fl1a = x_fl1a./S;
    
    cameraman.adcock = x_adcock;
    cameraman.fl1a = x_fl1a;
    cameraman.cs = x_cs;
    cameraman.l1a = x_l1a;
    
    cameraman.nbSamples = sum(sum(f_mask));
    cameraman.samplingRate = cameraman.nbSamples/m/n;
    cameraman.f_mask = f_mask;
    cameraman.masks = masks;
    
    cameraman.ssim.fl1a = ssim(x_fl1a,img_cam);
    cameraman.ssim.cs = ssim(x_cs,img_cam);
    cameraman.ssim.l1a = ssim(x_l1a,img_cam);
    cameraman.ssim.adcock = ssim(x_adcock,img_cam);
    
    cameraman.psnr.fl1a = psnr(x_fl1a,img_cam);
    cameraman.psnr.cs = psnr(x_cs,img_cam);
    cameraman.psnr.l1a = psnr(x_l1a,img_cam);
    cameraman.psnr.adcock = psnr(x_adcock,img_cam);
    
    cameraman.l2error.fl1a = norm(x_fl1a(:)-img_cam(:));
    cameraman.l2error.cs = norm(x_cs(:)-img_cam(:));
    cameraman.l2error.l1a = norm(x_l1a(:)-img_cam(:));
    cameraman.l2error.adcock = norm(x_adcock(:)-img_cam(:));
    
    figure
    subplot(2,4,1), imagesc(abs(x_cs)), colormap(gray(256)), title('CS Recovery')
    subplot(2,4,2), imagesc(abs(x_l1a)), colormap(gray(256)), title('L1 Analysis')
    subplot(2,4,3), imagesc(abs(x_fl1a)), colormap(gray(256)), title('Fused L1 analysis')
    subplot(2,4,4), imagesc(abs(x_adcock)), colormap(gray(256)), title('Central fusion')
    subplot(2,4,5), imagesc(abs(x_cs-img_cam)), colormap(gray(256)), title(['CS Error. SNR = ', num2str(psnr(x_cs,img_cam))])
    subplot(2,4,6), imagesc(abs(x_l1a-img_cam)), colormap(gray(256)), title(['L1 Error. SNR = ', num2str(psnr(x_l1a,img_cam))])
    subplot(2,4,7), imagesc(abs(x_fl1a-img_cam)), colormap(gray(256)), title(['Fused L1 Error. SNR = ', num2str(psnr(x_fl1a,img_cam))])
    subplot(2,4,8), imagesc(abs(x_adcock-img_cam)), colormap(gray(256)), title(['Central L1 Error. SNR = ', num2str(psnr(x_adcock,img_cam))])
    
end



%% Save everything

extensions = {'jpg', 'png', 'tif'};


vdsOutput_path = fullfile('.', 'results', 'CentralCompute', WAVEoptions.wname, [noiseName, num2str(noise_lvl)]);

if ~isfolder(vdsOutput_path)
    mkdir(vdsOutput_path);
end
 

for one_extension=1:length(extensions)
    % Start with the cameraman experiment -> camMan_dir_path
    % Save the sampling patterns
%     curRes = f_mask_I1; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_CScan_Radial.', extensions{one_extension}]))
%     curRes = f_mask_I2; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_CScan_Spiral.', extensions{one_extension}]))
%     curRes = f_mask_im1; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_RealBrain_Radial.', extensions{one_extension}]))
%     curRes = f_mask_im2; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['FourierPattern_RealBrain_Spiral.', extensions{one_extension}]))
%     
%     % Save the original images 
%     curRes = I; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['CScan_Original.', extensions{one_extension}]))
%     curRes = im; 
%     imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['RealBrain_Original.', extensions{one_extension}]))
    
    % Save the recoveries 
    % and the pointwise errors
    curRes = abs(x_adcock); 
    curErr = abs(x_adcock - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_Adcock.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_Adcock_error.', extensions{one_extension}]))
    curRes = abs(x_fl1a); 
    curErr = abs(x_fl1a - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_FL1A.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_FL1A_error.', extensions{one_extension}]))
    curRes = abs(x_l1a); 
    curErr = abs(x_l1a - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_L1A.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_L1A_error.', extensions{one_extension}]))
    curRes = abs(x_cs); 
    curErr = abs(x_cs - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_CS.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_CS_error.', extensions{one_extension}]))
    
    
end


save(fullfile(vdsOutput_path, 'allResults.m'))
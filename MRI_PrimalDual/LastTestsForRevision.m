% Tests used in the paper:
%
%
addpath('./sampling/') % Contains the various sampling schemes
addpath('./projections/') % Contains various filtering operations
addpath('./data/') % Some examples of data
addpath('./recovery/') % Some recovery procedures
addpath('./utils/') 

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
do_SL_FCS = true; % Run the tests for the phantom, with fused compressed sensing approach
do_SL_TV = true; % Run the tests on the phantom for Total Variation

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
    noise_lvl = 0.05;
    noiseName = 'Normal';
    % First for the multi-sensor case
    y_channel = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask, masks);
    % Then for the single sensor case
    y_rhs = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask);
    
    % Now we go ahead with the actual recovery
    x0 = [];
    measurement = @(f) c_fft_2d(f);
    measurementStar = @(f) c_ifft_2d(f);
%     x_adcock = SingleComputeFCS(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions, masks);
    x_cs = CS_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions);
    x_l1a = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    x_fl1a_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    % Compute the fused l1 analysis solution:
    x_fl1a = zeros(size(img_cam));
    for one_channel=1:nb_channels
        x_fl1a = x_fl1a + x_fl1a_all{one_channel};
    end
    x_fl1a = x_fl1a./S;
    
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
    
    cameraman.psnr.fl1a = psnr(x_fl1a,img_cam);
    cameraman.psnr.cs = psnr(x_cs,img_cam);
    cameraman.psnr.l1a = psnr(x_l1a,img_cam);
    
    cameraman.l2error.fl1a = norm(x_fl1a(:)-img_cam(:));
    cameraman.l2error.cs = norm(x_cs(:)-img_cam(:));
    cameraman.l2error.l1a = norm(x_l1a(:)-img_cam(:));
    
    figure
    subplot(2,3,1), imagesc(abs(x_cs)), colormap(gray(256)), title('CS Recovery')
    subplot(2,3,2), imagesc(abs(x_l1a)), colormap(gray(256)), title('L1 Analysis')
    subplot(2,3,3), imagesc(abs(x_fl1a)), colormap(gray(256)), title('Fused L1 analysis')
    subplot(2,3,4), imagesc(abs(x_cs-img_cam)), colormap(gray(256)), title('CS Error')
    subplot(2,3,5), imagesc(abs(x_l1a-img_cam)), colormap(gray(256)), title('L1 Error')
    subplot(2,3,6), imagesc(abs(x_fl1a-img_cam)), colormap(gray(256)), title('Fused L1 Error')
    
end


%% Phantom with fused approaches
if do_SL_FCS
    % Load the image
    N = 1024; 
    img_SL = phantom(N); % Shep-Logan CT Phantom
    [m,n] = size(img_SL);
    img_SL = img_SL/max(img_SL(:));
    
    % Generate sampling pattern
    how_spread = m^(1/2); % The bigger this parameter, the fewer sampling points will be taken
    f_mask = sample_Fourier_Gaussian(how_spread, m);
    
    % Generate sensor pre-filtering
    nb_channels = 8;
    masks = multi_channel_sines(m, nb_channels);
    S = zeros(size(img_SL));
    for one_channel=1:nb_channels
        S = S + masks{one_channel};
    end
    
    % Generate measurements
    noise_lvl = 0.05;
    noiseName = 'Normal';
    % First for the multi-sensor case
    y_channel = generate_noisy_Fourier_samples(img_SL, noise_lvl, noiseName, f_mask, masks);
    % Then for the single sensor case
    y_rhs = generate_noisy_Fourier_samples(img_SL, noise_lvl, noiseName, f_mask);
    
    % Now we go ahead with the actual recovery
    x0 = [];
    measurement = @(f) c_fft_2d(f);
    measurementStar = @(f) c_ifft_2d(f);
    x_l1a = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    x_fl1a_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    % Compute the fused l1 analysis solution:
    x_fl1a = zeros(size(img_SL));
    x_sos = zeros(size(img_SL));
    for one_channel=1:nb_channels
        x_fl1a = x_fl1a + x_fl1a_all{one_channel};
        x_sos = x_sos + x_fl1a_all{one_channel}.^2;
    end
    x_fl1a = x_fl1a./S;
    x_sos = sqrt(x_sos);
    
    SL_FCS.fl1a = x_fl1a;
    SL_FCS.sos = x_sos;
    SL_FCS.l1a = x_l1a;
    
    SL_FCS.nbSamples = sum(sum(f_mask));
    SL_FCS.samplingRate = SL_FCS.nbSamples/m/n;
    SL_FCS.f_mask = f_mask;
    SL_FCS.masks = masks;
    
    SL_FCS.ssim.fl1a = ssim(x_fl1a,img_SL);
    SL_FCS.ssim.sos = ssim(x_sos,img_SL);
    SL_FCS.ssim.l1a = ssim(x_l1a,img_SL);
    
    SL_FCS.psnr.fl1a = psnr(x_fl1a,img_SL);
    SL_FCS.psnr.sos = psnr(x_sos,img_SL);
    SL_FCS.psnr.l1a = psnr(x_l1a,img_SL);
    
    SL_FCS.l2error.fl1a = norm(x_fl1a(:)-img_SL(:));
    SL_FCS.l2error.sos = norm(x_sos(:)-img_SL(:));
    SL_FCS.l2error.l1a = norm(x_l1a(:)-img_SL(:));
    
    figure
    subplot(2,3,1), imagesc(abs(x_sos)), colormap(gray(256)), title('SOS Recovery')
    subplot(2,3,2), imagesc(abs(x_l1a)), colormap(gray(256)), title('L1 Analysis')
    subplot(2,3,3), imagesc(abs(x_fl1a)), colormap(gray(256)), title('Fused L1 analysis')
    subplot(2,3,4), imagesc(abs(x_sos-img_SL)), colormap(gray(256)), title('SOS Error')
    subplot(2,3,5), imagesc(abs(x_l1a-img_SL)), colormap(gray(256)), title('L1 Error')
    subplot(2,3,6), imagesc(abs(x_fl1a-img_SL)), colormap(gray(256)), title('Fused L1 Error')
    
end

%% Phantom with TV approaches

if do_SL_TV
    % Load the image
    N = 1024; 
    img_SL = phantom(N); % Shep-Logan CT Phantom
    [m,n] = size(img_SL);
    img_SL = img_SL/max(img_SL(:));
    
    % Generate sampling pattern
    how_spread = m^(1/2); % The bigger this parameter, the fewer sampling points will be taken
    f_mask = sample_Fourier_Gaussian(how_spread, m);
    
    % Generate sensor pre-filtering
    nb_channels = 8;
    % masks = multi_channel_sines(m, nb_channels);
    masks = multi_channel_lightning(n, nb_channels, [], n^2/4);
    S = zeros(size(img_SL));
    for one_channel=1:nb_channels
        S = S + masks{one_channel};
    end
    
    % Generate measurements
    noise_lvl = 0.05;
    noiseName = 'Normal';
    % First for the multi-sensor case
    y_channel = generate_noisy_Fourier_samples(img_SL, noise_lvl, noiseName, f_mask, masks);
    % Then for the single sensor case
    y_rhs = generate_noisy_Fourier_samples(img_SL, noise_lvl, noiseName, f_mask);
    
    % Now we go ahead with the actual recovery
    x0 = [];
    measurement = @(f) c_fft_2d(f);
    measurementStar = @(f) c_ifft_2d(f);
    x_fTV_all = FTV_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions);
    % x_TV = TV_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions);
    x_l1a = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    x_fl1a_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
    % Compute the fused l1 analysis solution:
    x_fl1a = zeros(size(img_SL));
    x_fTV = zeros(size(img_SL));
    x_sos = zeros(size(img_SL));
    for one_channel=1:nb_channels
        x_fl1a = x_fl1a + x_fl1a_all{one_channel};
        x_fTV = x_fTV + x_fTV_all{one_channel};
        x_sos = x_sos + x_fl1a_all{one_channel}.^2;
    end
    
    x_sos = sqrt(x_sos);
    x_fl1a = x_fl1a./S;
    x_fTV = x_fTV./S;
    
    SL_TV.fl1a = x_fl1a;
    SL_TV.fTV = x_fTV;
    SL_TV.l1a = x_l1a;
    SL_TV.sos = x_sos;
    
    SL_TV.nbSamples = sum(sum(f_mask));
    SL_TV.samplingRate = SL_TV.nbSamples/m/n;
    SL_TV.f_mask = f_mask;
    SL_TV.masks = masks;
    
    SL_TV.ssim.fl1a = ssim(x_fl1a,img_SL);
    SL_TV.ssim.fTV = ssim(x_fTV,img_SL);
    SL_TV.ssim.l1a = ssim(x_l1a,img_SL);
    SL_TV.ssim.sos = ssim(x_sos,img_SL);
    
    SL_TV.psnr.fl1a = psnr(x_fl1a,img_SL);
    SL_TV.psnr.fTV = psnr(x_fTV,img_SL);
    SL_TV.psnr.l1a = psnr(x_l1a,img_SL);
    SL_TV.psnr.sos = psnr(x_sos,img_SL);
    
    
    SL_TV.l2error.fl1a = norm(x_fl1a(:)-img_SL(:));
    SL_TV.l2error.fTV = norm(x_fTV(:)-img_SL(:));
    SL_TV.l2error.l1a = norm(x_l1a(:)-img_SL(:));
    SL_TV.l2error.sos = norm(x_sos(:) - img_SL(:));
    
    figure
    subplot(2,4,1), imagesc(abs(x_fTV)), colormap(gray(256)), title('Fused TV Recovery')
    subplot(2,4,2), imagesc(abs(x_l1a)), colormap(gray(256)), title('L1 Analysis')
    subplot(2,4,3), imagesc(abs(x_fl1a)), colormap(gray(256)), title('Fused L1 analysis')
    subplot(2,4,4), imagesc(abs(x_sos)), colormap(gray(256)), title('SOS Recovery')
    subplot(2,4,5), imagesc(abs(x_fTV-img_SL)), colormap(gray(256)), title('fTV Error')
    subplot(2,4,6), imagesc(abs(x_l1a-img_SL)), colormap(gray(256)), title('L1 Error')
    subplot(2,4,7), imagesc(abs(x_fl1a-img_SL)), colormap(gray(256)), title('Fused L1 Error')
    subplot(2,4,8), imagesc(abs(x_sos-img_SL)), colormap(gray(256)), title('SOS Error')
    
    
end


%% Saving everything
extensions = {'jpg', 'png', 'tif'};

camMan_dir_path = fullfile('.','results', 'originalPaper', 'cameraman', ['samples', num2str(num2str(cameraman.nbSamples))], [noiseName, num2str(noise_lvl)] );
SL_FCS_dir_path = fullfile('.','results', 'originalPaper', 'phantom_sines', ['samples', num2str(num2str(SL_FCS.nbSamples))], [noiseName, num2str(noise_lvl)] );
SL_TV_dir_path = fullfile('.','results', 'originalPaper', 'phantom_beams', ['samples', num2str(num2str(SL_TV.nbSamples))], [noiseName, num2str(noise_lvl)] );

if ~isfolder(camMan_dir_path)
    mkdir(camMan_dir_path)
end
if ~isfolder(SL_FCS_dir_path)
    mkdir(SL_FCS_dir_path)
end 
if ~isfolder(SL_TV_dir_path)
    mkdir(SL_TV_dir_path)
end 

for one_extension=1:length(extensions)
    % Start with the cameraman experiment -> camMan_dir_path
    % CS recovery
    curRes = abs(cameraman.cs); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(cameraman.cs - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(camMan_dir_path, ['camera_CS_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(camMan_dir_path, ['camera_CS_error.', extensions{one_extension}]))
    % Fused l1 analysis
    curRes = abs(cameraman.fl1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(cameraman.fl1a - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(camMan_dir_path, ['camera_fL1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(camMan_dir_path, ['camera_fL1A_error.', extensions{one_extension}]))
    % Plain l1 analysis
    curRes = abs(cameraman.l1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(cameraman.l1a - img_cam);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(camMan_dir_path, ['camera_L1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(camMan_dir_path, ['camera_L1A_error.', extensions{one_extension}]))
    
    % We continue with the phantom with sines -> SL_FCS_dir_path
    % Sum of Squared recovery
    curRes = abs(SL_FCS.sos); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_FCS.sos - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_SOS_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_SOS_error.', extensions{one_extension}]))
    % Fused l1 analysis
    curRes = abs(SL_FCS.fl1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_FCS.fl1a - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_fL1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_fL1A_error.', extensions{one_extension}]))
    % Plain l1 analysis
    curRes = abs(SL_FCS.l1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_FCS.l1a - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_L1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_FCS_dir_path, ['SL_FCS_L1A_error.', extensions{one_extension}]))
    
    % And finally the phantom with spherical beams -> SL_TV_dir_path
    % Sum of Squared recovery
    curRes = abs(SL_TV.sos); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_TV.sos - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_TV_dir_path, ['SL_TV_SOS_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_TV_dir_path, ['SL_TV_SOS_error.', extensions{one_extension}]))
    % Fused l1 analysis
    curRes = abs(SL_TV.fl1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_TV.fl1a - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_TV_dir_path, ['SL_TV_fL1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_TV_dir_path, ['SL_TV_fL1A_error.', extensions{one_extension}]))
    % Plain l1 analysis
    curRes = abs(SL_TV.l1a); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_TV.l1a - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_TV_dir_path, ['SL_TV_L1A_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_TV_dir_path, ['SL_TV_L1A_error.', extensions{one_extension}]))
    % Fused TV minimization
    curRes = abs(SL_TV.fTV); % Current result, the one being saved! Abs is to make sure that the true 0 is black
    curErr = abs(SL_TV.fTV - img_SL);
    imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(SL_TV_dir_path, ['SL_TV_fTV_results.', extensions{one_extension}]))
    imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(SL_TV_dir_path, ['SL_TV_fTV_error.', extensions{one_extension}]))
    
   
end

save(fullfile('.','results', 'originalPaper', 'allResults.m'))
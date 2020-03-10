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

preMesurementNoises = 0:0.001:0.026;

%% Cameraman image
if do_cameraman
    % Load the image
    img_cam = double(imread('cameraman.tif'));
    [m,n] = size(img_cam);
    img_cam = img_cam/max(img_cam(:));
    noiseLessImg = img_cam;
    
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
    
    noise_lvl = 0.01;
        noiseName = 'Normal';
    
    for oneNoise_idx=1:length(preMesurementNoises)
        current_noise = preMesurementNoises(oneNoise_idx);
        img_cam = imnoise(noiseLessImg, 'salt & pepper', current_noise);
        
        % Generate measurements
        % First for the multi-sensor case
        y_channel = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask, masks);
        % Then for the single sensor case
        y_rhs = generate_noisy_Fourier_samples(img_cam, noise_lvl, noiseName, f_mask);
        
        % Now we go ahead with the actual recovery
        x0 = [];
        measurement = @(f) c_fft_2d(f);
        measurementStar = @(f) c_ifft_2d(f);
        x_cs = CS_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions);
        x_l1a = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
        x_fl1a_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions);
        x_fcs_all = FL1A_recovery(y_channel,f_mask, x0, measurement, measurementStar, PDoptions, []);
        % Compute the fused l1 analysis solution:
        x_fl1a = zeros(size(img_cam));
        x_fcs = zeros(size(img_cam));
        for one_channel=1:nb_channels
            x_fl1a = x_fl1a + x_fl1a_all{one_channel};
            x_fcs = x_fcs + x_fcs_all{one_channel};
        end
        x_fl1a = x_fl1a./S;
        x_fcs = x_fcs./S;
        
        cameraman.fcs = x_fcs;
        cameraman.fl1a = x_fl1a;
        cameraman.cs = x_cs;
        cameraman.l1a = x_l1a;
        
        cameraman.nbSamples = sum(sum(f_mask));
        cameraman.samplingRate = cameraman.nbSamples/m/n;
        cameraman.f_mask = f_mask;
        cameraman.masks = masks;
        
        cameraman.ssim.fl1a = ssim(x_fl1a,noiseLessImg);
        cameraman.ssim.cs = ssim(x_cs,noiseLessImg);
        cameraman.ssim.l1a = ssim(x_l1a,noiseLessImg);
        cameraman.ssim.fcs = ssim(x_fcs,noiseLessImg);
        
        cameraman.psnr.fl1a = psnr(x_fl1a,noiseLessImg);
        cameraman.psnr.cs = psnr(x_cs,noiseLessImg);
        cameraman.psnr.l1a = psnr(x_l1a,noiseLessImg);
        cameraman.psnr.fcs = psnr(x_fcs,noiseLessImg);
        
        cameraman.l2error.fl1a = norm(x_fl1a(:)-noiseLessImg(:));
        cameraman.l2error.cs = norm(x_cs(:)-noiseLessImg(:));
        cameraman.l2error.l1a = norm(x_l1a(:)-noiseLessImg(:));
        cameraman.l2error.fcs = norm(x_fcs(:)-noiseLessImg(:));
        
        figure
        subplot(2,4,1), imagesc(abs(x_cs)), colormap(gray(256)), title('CS Recovery')
        subplot(2,4,2), imagesc(abs(x_l1a)), colormap(gray(256)), title('L1 Analysis')
        subplot(2,4,3), imagesc(abs(x_fl1a)), colormap(gray(256)), title('Fused L1 analysis')
        subplot(2,4,4), imagesc(abs(x_fcs)), colormap(gray(256)), title('Fused CS')
        subplot(2,4,5), imagesc(abs(x_cs-noiseLessImg)), colormap(gray(256)), title(['CS Error. SNR = ', num2str(psnr(x_cs,noiseLessImg))])
        subplot(2,4,6), imagesc(abs(x_l1a-noiseLessImg)), colormap(gray(256)), title(['L1 Error. SNR = ', num2str(psnr(x_l1a,noiseLessImg))])
        subplot(2,4,7), imagesc(abs(x_fl1a-noiseLessImg)), colormap(gray(256)), title(['Fused L1 Error. SNR = ', num2str(psnr(x_fl1a,noiseLessImg))])
        subplot(2,4,8), imagesc(abs(x_fcs-noiseLessImg)), colormap(gray(256)), title(['Fused CS Error. SNR = ', num2str(psnr(x_fcs,noiseLessImg))])
        
        
        
        
        
        %% Save everything
        
        extensions = {'jpg', 'png', 'tif'};
        
        
        vdsOutput_path = fullfile('.', 'results', 'NoiseBehaviour', WAVEoptions.wname, ['CurrNoise', num2str(current_noise)]);
        
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
            curRes = abs(x_fcs);
            curErr = abs(x_fcs - noiseLessImg);
            imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_FCS.', extensions{one_extension}]))
            imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_FCS_error.', extensions{one_extension}]))
            curRes = abs(x_fl1a);
            curErr = abs(x_fl1a - noiseLessImg);
            imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_FL1A.', extensions{one_extension}]))
            imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_FL1A_error.', extensions{one_extension}]))
            curRes = abs(x_l1a);
            curErr = abs(x_l1a - noiseLessImg);
            imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_L1A.', extensions{one_extension}]))
            imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_L1A_error.', extensions{one_extension}]))
            curRes = abs(x_cs);
            curErr = abs(x_cs - noiseLessImg);
            imwrite((curRes-min(curRes(:)))/(max(curRes(:))-min(curRes(:))),fullfile(vdsOutput_path, ['Cameraman_CS.', extensions{one_extension}]))
            imwrite((curErr-min(curErr(:)))/(max(curErr(:))-min(curErr(:))),fullfile(vdsOutput_path, ['Cameraman_CS_error.', extensions{one_extension}]))
            
            curImg = img_cam;
            imwrite((curImg-min(curImg(:)))/(max(curImg(:))-min(curImg(:))),fullfile(vdsOutput_path, ['SnPNoisyImage.', extensions{one_extension}]))
            
        end
        
        
        save(fullfile(vdsOutput_path, 'allResults.m'))
        
    end
end
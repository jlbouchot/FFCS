function [l2error, snrValues, ssimValues, dir_path] = saveResults(x_hats, img, experiment, options)


%% Check do_XX
nb_doers = 9;
doers = {'do_cs', 'do_sos','do_fcs', 'do_ftv', 'do_l1s', 'do_l1a', 'do_tv', 'do_fl1s', 'do_fl1a'};
for one_doer_idx=1:nb_doers
    cur_doer = doers{one_doer_idx};
    if ~isfield(options,cur_doer) % makes sure we haven't forgotten to set a certain field
        options = setfield(options, cur_doer, false);
    end
end

%% Some preparation
[m,n] = size(x_hats.fcs);
wave_name = options.wname; 
sampling_ratio = sum(options.f_mask(:))/m/n;
nbSamples = sum(options.f_mask(:));

x_cs = x_hats.cs; 
x_fcs = x_hats.fcs;
x_TV = x_hats.tv; 
x_ftv = x_hats.ftv;
x_l1s = x_hats.l1s;
x_fl1s = x_hats.fl1s;
x_l1a = x_hats.l1a;
x_fl1a = x_hats.fl1a;
x_sos = x_hats.sos;

%% Compute pointwise errors for all situations
err_cs = abs(img-x_cs);
l2error.cs = norm(err_cs);

err_fcs = abs(img-x_fcs);
l2error.fcs = norm(err_fcs);

err_tv = abs(img-x_TV);
l2error.tv = norm(err_tv);

err_ftv = abs(img-x_ftv);
l2error.ftv = norm(err_ftv);

err_l1a = abs(img-x_l1a);
l2error.l1a = norm(err_l1a);

err_fl1a = abs(img-x_fl1a);
l2error.fl1a = norm(err_fl1a);

err_l1s = abs(img-x_l1s);
l2error.l1s = norm(err_l1s);

err_fl1s = abs(img-x_fl1s);
l2error.fl1s = norm(err_fl1s);

err_sos = abs(img - x_sos);
l2error.sos = norm(err_sos)

%% Save results and images
output_folder = fullfile('.', 'results')
experiment
experiment_name = fullfile(experiment, ['N', num2str(n)], ['wave', wave_name], ['nbSamples', num2str(nbSamples)], [options.noiseName, num2str(options.noise_lvl)]);

dir_path = fullfile(output_folder, experiment_name)

if isfolder(dir_path) == false
    mkdir(dir_path)
end

extensions = {'jpg', 'png', 'tif'};

for one_extension=1:length(extensions)
    imwrite((x_fcs-min(x_fcs(:)))/(max(x_fcs(:))-min(x_fcs(:))),fullfile(dir_path, ['FCS_results.', extensions{one_extension}]))
    imwrite((x_ftv-min(x_ftv(:)))/(max(x_ftv(:))-min(x_ftv(:))),fullfile(dir_path, ['FCS_TV.', extensions{one_extension}]))
    imwrite((x_fl1a-min(x_fl1a(:)))/(max(x_fl1a(:))-min(x_fl1a(:))),fullfile(dir_path, ['Fused-L1-Analysis.', extensions{one_extension}]))
    imwrite((x_l1a-min(x_l1a(:)))/(max(x_l1a(:))-min(x_l1a(:))),fullfile(dir_path, ['L1-Analysis.', extensions{one_extension}]))
    imwrite((x_fl1s-min(x_fl1s(:)))/(max(x_fl1s(:))-min(x_fl1s(:))),fullfile(dir_path, ['Fused-L1-Synthesis.', extensions{one_extension}]))
    imwrite((x_l1s-min(x_l1s(:)))/(max(x_l1s(:))-min(x_l1s(:))),fullfile(dir_path, ['L1-Synthesis.', extensions{one_extension}]))
    imwrite((x_TV-min(x_TV(:)))/(max(x_TV(:))-min(x_TV(:))),fullfile(dir_path, ['TV-Min.', extensions{one_extension}]))
    imwrite((x_sos-min(x_sos(:)))/(max(x_sos(:))-min(x_sos(:))),fullfile(dir_path, ['Sum-of-Squares.', extensions{one_extension}]))
    imwrite((x_cs-min(x_cs(:)))/(max(x_cs(:))-min(x_cs(:))),fullfile(dir_path, ['Traditional_CS_results.', extensions{one_extension}]))
    
    % Pointwise difference for CS
    imwrite((err_cs - min(err_cs(:))) / (max(err_cs(:)) - min(err_cs(:)) ), fullfile(dir_path, ['Traditional_CS_error.', extensions{one_extension}])) 
    % Pointwise difference for TV
    imwrite((err_tv - min(err_tv(:))) / (max(err_tv(:)) - min(err_tv(:)) ), fullfile(dir_path, ['TV-min_error.', extensions{one_extension}])) 
    % Pointwise difference for l1 analysis
    imwrite((err_l1a - min(err_l1a(:))) / (max(err_l1a(:)) - min(err_l1a(:)) ), fullfile(dir_path, ['L1-Analysis_error.', extensions{one_extension}])) 
    % Pointwise difference for fused l1 analysis
    imwrite((err_fl1a - min(err_fl1a(:))) / (max(err_fl1a(:)) - min(err_fl1a(:)) ), fullfile(dir_path, ['Fused-L1-Analysis_error.', extensions{one_extension}])) 
    % Pointwise difference for l1 analysis
    imwrite((err_l1s - min(err_l1s(:))) / (max(err_l1s(:)) - min(err_l1s(:)) ), fullfile(dir_path, ['L1-Synthesis_error.', extensions{one_extension}])) 
    % Pointwise difference for fused l1 analysis
    imwrite((err_fl1s - min(err_fl1s(:))) / (max(err_fl1s(:)) - min(err_fl1s(:)) ), fullfile(dir_path, ['Fused-L1-Synthesis_error.', extensions{one_extension}])) 
    % Pointwise difference for fcs
    imwrite((err_fcs - min(err_fcs(:))) / (max(err_fcs(:)) - min(err_fcs(:)) ), fullfile(dir_path, ['FCS_error.', extensions{one_extension}])) 
    % Pointwise difference for fcs with TV
    imwrite((err_ftv - min(err_ftv(:))) / (max(err_ftv(:)) - min(err_ftv(:)) ), fullfile(dir_path, ['FCS_TV_error.', extensions{one_extension}])) 
    % Pointwise difference for sos
    imwrite((err_sos - min(err_sos(:))) / (max(err_sos(:)) - min(err_sos(:)) ), fullfile(dir_path, ['Sum-of-Squares_error.', extensions{one_extension}])) 
end

%% Compute SNRs
snrValues.cs = psnr(x_cs,img);
snrValues.l1a = psnr(x_l1a,img);
snrValues.fl1a = psnr(x_fl1a,img);
snrValues.l1s = psnr(x_l1s,img);
snrValues.fl1s = psnr(x_fl1s,img);
snrValues.tv = psnr(x_TV,img);
snrValues.sos = psnr(x_sos,img);
snrValues.fcs = psnr(x_fcs,img);
snrValues.fcs_TV = psnr(x_ftv,img)

%% Compute SSIMs
ssimValues.cs = ssim(x_cs,img);
ssimValues.l1a = ssim(x_l1a,img);
ssimValues.fl1a = ssim(x_fl1a,img);
ssimValues.l1s = ssim(x_l1s,img);
ssimValues.fl1s = ssim(x_fl1s,img);
ssimValues.tv = ssim(x_TV,img);
ssimValues.sos = ssim(x_sos,img);
ssimValues.fcs = ssim(x_fcs,img);
ssimValues.fcs_TV = ssim(x_ftv,img)

%% Save .eps files
% Have to deal differently for .eps files

imagesc(x_fcs);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'FCS_results.eps') ,'-deps2','-r300');  %# Print the figure


imagesc(x_cs);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'CS_results.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'CS_results.eps'],'-deps2','-r300');  %# Print the figure



imagesc(x_ftv);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'FCS_TV.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'FCS_TV.eps'],'-deps2','-r300');  %# Print the figure

imagesc(x_l1s);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'L1-Synthesis.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'L1-Synthesis.eps'],'-deps2','-r300');  %# Print the figure

imagesc(x_fl1s);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'Fused-L1-Synthesis.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'Fused-L1-Synthesis.eps'],'-deps2','-r300');  %# Print the figure

imagesc(x_l1a);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'L1-Analysis.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'L1-Analysis.eps'],'-deps2','-r300');  %# Print the figure

imagesc(x_fl1a);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'Fused-L1-Analysis.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, 'Fused-L1-Analysis.eps'],'-deps2','-r300');  %# Print the figure

imagesc(x_sos);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'Sum-of-Squares.eps') ,'-deps2','-r300');  %# Print the figure
% print(gcf,[dir_path, ''],'-deps2','-r300');  %# Print the figure

imagesc(x_TV);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,fullfile(dir_path, 'TV-min.eps') ,'-deps2','-r300');  %# Print the figure

close all

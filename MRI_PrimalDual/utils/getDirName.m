function [outputFolder] = getDirName(experiment, m, n, options)

wave_name = options.wname; 
sampling_ratio = sum(options.f_mask(:))/m/n;
nbSamples = sum(options.f_mask(:));

output_folder = fullfile('.', 'results');

experiment_name = fullfile(experiment, ['N', num2str(n)], ['wave', wave_name], ['nbSamples', num2str(nbSamples)], ['noise', num2str(options.noise_lvl)]);

outputFolder = fullfile(output_folder, experiment_name);
end
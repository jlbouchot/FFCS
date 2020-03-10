% Using our fused compressed sensing approach to recover noisy MRI 

%% Bit of cleaning
clear all
close all
clc

%% Add the important files 
addpath('./sampling/') % Contains the various sampling schemes
addpath('./projections/') % Contains various filtering operations
addpath('./data/') % Some examples of data

%% Experiment 1: TV minimization, with subsampled Fourier measurements
% We consider here :
% * Artificial image phantom()
% * Size N = 1024 (can be changed easily!)
% * Subsampled Fourier samples, according to a Gaussian distribution
% * The minimization is performed with a self implemented primal dual
% algorithm
% * 
% 
% /!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\/!\
% ------ WARNINGS
% * The way the sampling ratio is done is dependent on the size. This could
% and should be done better. 
% * For the proximal operator needed in the primal dual algorithm, we have
% used a first order approximation of a Neumann series to approximate the
% inverse of I + tau F*F. With tau small enough, this is already pretty
% good. 
% 
% The results are compared against
% * Basic compressed sensing
% * TV minimization
% * Simple sum of squares combination with Daubechies wavelets

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters used everywhere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In l1 analysis compressed sensing, we minize the l1 norm of the analysis
% coefficients. Astar is the synthesis operator. 
wave_name = 'db4';
nb_level = 4;

% Generate the image
artificial = 1; % Do we want an artificial model?
if artificial == 1
    m = 256; n = m;
    img = phantom(n); % Shep-Logan CT Phantom
    
%     img = img/max(img(:));
    
    experiment = 'artificial'; % used for automatic saving
elseif artificial == 2
    load brain;
    [m,n] = size(im);
    img = real(im);
    clear im;
    
    experiment = 'real'; % used for automatic saving
else 
    img = double(imread('cameraman.tif'));
    [m,n] = size(img);    
    img = img/max(img(:));
    
    experiment = 'camera';
end




% Parameters for the Primal dual algorithm
sigma = 1; 
tau=.95/16/sigma; 
theta=1;
nb_max_iter = 400;

% Add a bit of noise? 
noise_lvl = 0.0; % Noise in the measurement process
img_noise = 0.0; % Noise in the input image

original_img = img; 
img = img + img_noise*randn(m,n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the subsampling pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose from "radial", "line", "uniform", "gaussian", "horizontal", "vertical"
sampling_type = "gaussian";


switch lower(sampling_type)
    case 'radial'
        % Lines through the center of the Fourier domain
        nb_lines = ceil(m/12); % This corresponds to roughly 8%
        f_mask = sample_Fourier_angles(nb_lines, m);
    case "lines"
        % Lines parallel to the canonical directions
        nb_line_per_dim = ceil(m/8); % This corresponds to roughly 12%
        f_mask = sample_Fourier_lines(nb_line_per_dim, m);
    case "uniform"
        % Uniform Fourier sampling
        rejection_threshold = 0.1;
        f_mask = sample_Fourier_uniform(rejection_threshold, m);
    case "horizontal"
        nb_lines = m/5;
        f_mask = sample_Fourier_oneD(nb_lines, m, false);
    case "vertical"
        nb_lines = m/7;
        f_mask = sample_Fourier_oneD(nb_lines, m, true);
    case "gaussian"
        % Gaussian centered sampling
        how_spread = m^(1/2); % The bigger this parameter, the fewer sampling points will be taken
        f_mask = sample_Fourier_Gaussian(how_spread, m);
    case "VDSCentre"
        % Variable density sampling with deterministic sampling in the
        % centre. The outside will be considered as radial sampling
        
    otherwise
        disp("FCS-MRI Error: Did not understand your sampling type. Will use Gaussian");
        % Gaussian centered sampling
        how_spread = m^(1/2); % The bigger this parameter, the fewer sampling points will be taken
        f_mask = sample_Fourier_Gaussian(how_spread, m);
end


sampling_ratio = num2str(sum(sum(f_mask))/numel(f_mask)*100);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-channel options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % Random Gaussian-like illumination things
% % nb_channels = 10; % This will also corresponds to the number of stations. 
% % masks = multi_channel_lightning(n, nb_channels, [], n^2/4);


% % Using Gaussian-like illumination things
% centers = [floor(n/8), floor(n/8); floor(n/8), floor(7*n/8); floor(7*n/8), floor(n/8); floor(7*n/8), floor(7*n/8); floor(n/2), floor(n/2)];
% nb_channels = size(centers,1); % This will also corresponds to the number of stations. 
% masks = multi_channel_lightning(n, nb_channels, centers, n^2);


% % % % Using random canonical projections
% % % nb_channels = 15; 
% % % rank = n 
% % % [masks, S] = random_canonical_proj(nb_proj, N, rank)

% % % Splitting of the original domain
% % masks = multi_channel_hard_assignments(m);
% % nb_channels = length(masks);
% % disp(num2str(nb_channels))

% Generate a field of sine-weighted illuminations
nb_channels = 8;
masks = multi_channel_sines(m, nb_channels);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Traditional compressed sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using traditional compressed sensing... ")

disp(['Recovery from ', sampling_ratio, ' % measurements'])

% Get the measurements
img_f = c_fft_2d(img) + noise_lvl*randn(m,n);
y_rhs = f_mask.*img_f;


% In traditional compressed sensing, there are no operators inside the
% first norm. Hence, A and Astar are simply the identity. 
A = @(f) f;
Astar = @(f) f; 

% Proximal operator for the l1 norm
prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
% Proximal operator for the l2 data fidelity term
prox_G = @(z) (real(c_ifft_2d(y_rhs)) + z) - ... 
    tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_rhs) + z)));

% Call the primal dual algorithm
disp('          ... primal dual algorithm running ... ')
cs_start = tic;
x_cs = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
cs_stop = toc(cs_start)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Traditional compressed sensing with Douglas-Rachford
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using traditional compressed sensing and DRS... ")

disp(['Recovery from ', sampling_ratio, ' % measurements'])

% Get the measurements
img_f = c_fft_2d(img) + noise_lvl*randn(m,n);
y_rhs = f_mask.*img_f;


% Here A and Astar are seen as measurements and its adjoint
A = @(f) f_mask.*c_fft_2d(f);
Astar = @(f) c_ifft_2d(f); 

TV = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
TVstar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));

% Proximal operator for the l1 norm
prox_F = @(tau, z) sign(z).*(abs(z)-tau).*(abs(z)>=tau) ;
% Proximal operator for the l2 data fidelity term
prox_G = @(z, tau) z +real(Astar(y_rhs-A(z)));

% Call the primal dual algorithm
disp('          ... Douglas-Rachford Splitting in progress ... ')
cs_start = tic;
[x_DRS, l1values] = DouglasRachfordMin(zeros(m,n), prox_F, prox_G, 10*nb_max_iter, 1, 0.01, 50);
% PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
cs_stop = toc(cs_start)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TV minimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery by total variation minimization... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])

% Get the measurements
img_f = c_fft_2d(img) + noise_lvl*randn(m,n);
y_rhs = f_mask.*img_f;


A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));

abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

prox_G = @(z) (real(c_ifft_2d(y_rhs)) + z) - ... 
    tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_rhs) + z)));
prox_G = @(z) z + real(c_ifft_2d(y_rhs - f_mask.*c_fft_2d(z))); 

% Call the primal dual algorithm
disp('          ... primal dual algorithm running ... ')
nb_max_iter = 200;
TV_start = tic;
x_TV = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
TV_stop = toc(TV_start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 analysis compressed sensing (Daubechies 4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using l1 analysis... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])

% Get the measurements
img_f = c_fft_2d(img) + noise_lvl*randn(m,n);
y_rhs = f_mask.*img_f;

% Get the structure first, this will not change through time
[~,s]=wavedec2(img,nb_level,wave_name);

% Analysis operator
A = @(f) wavedec2(f,nb_level,wave_name);
% Synthesis operator
Astar = @(f) waverec2(f,s,wave_name);

prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

prox_G = @(z) (real(c_ifft_2d(y_rhs)) + z) - ... 
    tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_rhs) + z)));

% Call the primal dual algorithm
disp('          ... primal dual algorithm running ... ')
l1_start = tic;
x_l1 = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
l1_stop = toc(l1_start)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1 analysis compressed sensing -- sum of squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using l1 analysis SoS ... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])


disp([num2str(nb_channels), ' channels will be used. '])

% Get the measurements
for one_channel=1:nb_channels
    img_channels{one_channel} = c_fft_2d(masks{one_channel}.*img) + noise_lvl*randn(m,n);
    y_channel{one_channel} = f_mask.*img_channels{one_channel};
end


% Get the structure first, this will not change through time
[~,s]=wavedec2(img_channels{1},nb_level,wave_name);
% Analysis operator
A = @(f) wavedec2(f,nb_level,wave_name);
% Synthesis operator
Astar = @(f) waverec2(f,s,wave_name);

% Initialize the end result with zeros
x_sos = zeros(m,n);

sos_stops = zeros(nb_channels+1,1);
for one_channel=1:nb_channels
    sos_start = tic;
    % Proximal operator for the l1 norm
    prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
    % Proximal operator for the l2 data fidelity term
    %prox_G = @(z) (real(c_ifft_2d(y_channel{one_channel})) + z) - ...
    %    tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_channel{one_channel}) + z)));
    prox_G = @(z) z + real(c_ifft_2d(y_channel{one_channel} - f_mask.*c_fft_2d(z))); % + (real(c_ifft_2d(y_channel{one_channel})) + z) - ...
       
    
    % Call the primal dual algorithm
    disp(['          ... primal dual algorithm running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
    x_local{one_channel} = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
    x_sos = x_sos + x_local{one_channel}.^2;
    sos_stops(one_channel) = toc(sos_start);
end

sos_start = tic;
x_sos = sqrt(x_sos);
sos_stops(nb_channels+1) = toc(sos_start);

sos_stop = sum(sos_stops)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fused L1 analysis compressed sensing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using l1 analysis FCS ... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])




disp([num2str(nb_channels), ' channels will be used. '])

% Get the measurements
for one_channel=1:nb_channels
    img_channels{one_channel} = c_fft_2d(masks{one_channel}.*img) + noise_lvl*randn(m,n);
    y_channel{one_channel} = f_mask.*img_channels{one_channel};
end



% Get the structure first, this will not change through time
[~,s]=wavedec2(img_channels{1},nb_level,wave_name);
% Analysis operator
A = @(f) wavedec2(f,nb_level,wave_name);
% Synthesis operator
Astar = @(f) waverec2(f,s,wave_name);

% Initialize the end result with zeros
x_fcs = zeros(m,n);

fcs_stops = zeros(nb_channels+1,1);
for one_channel=1:nb_channels
    fcs_start = tic;
    % Proximal operator for the l1 norm
    prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
    % Proximal operator for the l2 data fidelity term
%     prox_G = @(z) (real(c_ifft_2d(y_channel{one_channel})) + z) - ...
%         tau*real(c_ifft_2d(f_mask.*c_fft_2d(c_ifft_2d(y_channel{one_channel}) + z)));
    prox_G = @(z) z + real(c_ifft_2d(y_channel{one_channel} - f_mask.*c_fft_2d(z))); 
    
    % Call the primal dual algorithm
    disp(['          ... primal dual algorithm running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
    x_local{one_channel} = PrimalDualMin(zeros(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
    x_fcs = x_fcs + x_local{one_channel};
    fcs_stops(one_channel) = toc(fcs_start);
end

fcs_start = tic;
S = zeros(m,n);
for i=1:nb_channels
    S = S + masks{i};
end
x_fcs = x_fcs./S;
fcs_stops(nb_channels+1) = toc(fcs_start); 
fcs_stop = sum(fcs_stops)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused TV mminimization 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using TV-FCS ... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])


disp([num2str(nb_channels), ' channels will be used. '])

% Get the measurements
for one_channel=1:nb_channels
    img_channels{one_channel} = c_fft_2d(masks{one_channel}.*img) + noise_lvl*randn(m,n);
    y_channel{one_channel} = f_mask.*img_channels{one_channel};
end

A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));

abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

% Initialize the end result with zeros
x_fcs_TV = zeros(m,n);

fcs_TV_stops = zeros(nb_channels+1,1);
for one_channel=1:nb_channels
    fcs_TV_start = tic;
    % Proximal operator for the l2 data fidelity term
    prox_G = @(z) (real(c_ifft_2d(y_channel{one_channel})) + z) - ...
        tau*real(c_ifft_2d(f_mask.*c_fft_2d(masks{one_channel}.*c_ifft_2d(y_channel{one_channel}) + z)));
    
    % Call the primal dual algorithm
    disp(['          ... primal dual algorithm running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
    x_local{one_channel} = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
    x_fcs_TV = x_fcs_TV + x_local{one_channel};
    fcs_TV_stops(one_channel) = toc(fcs_TV_start);
end

fcs_TV_start = tic;
S = zeros(m,n);
for i=1:nb_channels
    S = S + masks{i};
end
x_fcs_TV = x_fcs_TV./S;
fcs_TV_stops(nb_channels+1) = toc(fcs_TV_start);
fcs_TV_stop = sum(fcs_TV_stops)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fused TV mminimization with other G function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp("Computing the recovery using TV-FCS number 2!... ")
disp(['Recovery from ', sampling_ratio, ' % measurements'])


disp([num2str(nb_channels), ' channels will be used. '])

% Get the measurements
for one_channel=1:nb_channels
    img_channels{one_channel} = c_fft_2d(masks{one_channel}.*img) + noise_lvl*randn(m,n);
    y_channel{one_channel} = f_mask.*img_channels{one_channel};
end

A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));

abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

% Initialize the end result with zeros
x_fcs_TV2 = zeros(m,n);

for one_channel=1:nb_channels
    % Proximal operator for the l2 data fidelity term
%     prox_G = @(z) (real(c_ifft_2d(y_channel{one_channel})) + z) - ...
%         tau*real(c_ifft_2d(f_mask.*c_fft_2d(masks{one_channel}.*c_ifft_2d(y_channel{one_channel}) + z)));
    prox_G = @(z) z + real(c_ifft_2d(y_channel{one_channel} - z)); 
    
    % Call the primal dual algorithm
    disp(['          ... primal dual algorithm running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
    x_local{one_channel} = PrimalDualMin(ones(m,n), prox_Fstar, prox_G, A, Astar, nb_max_iter, theta, sigma, tau);
    x_fcs_TV2 = x_fcs_TV2 + x_local{one_channel};
end

S = zeros(m,n);
for i=1:nb_channels
    S = S + masks{i};
end
x_fcs_TV2 = x_fcs_TV2./S;



%% Some display
figure
colormap(gray(256))
imagesc([(x_cs-min(x_cs(:)))/(max(x_cs(:))-min(x_cs(:))), (x_l1-min(x_l1(:)))/(max(x_l1(:))-min(x_l1(:))), (x_TV-min(x_TV(:)))/(max(x_TV(:))-min(x_TV(:))); (x_sos-min(x_sos(:)))/(max(x_sos(:))-min(x_sos(:))), (x_fcs-min(x_fcs(:)))/(max(x_fcs(:))-min(x_fcs(:))), (x_fcs_TV-min(x_fcs_TV(:)))/(max(x_fcs_TV(:))-min(x_fcs_TV(:)))])


figure
colormap(gray(256))
imagesc([masks{1}.*img,masks{2}.*img,masks{3}.*img;masks{4}.*img,masks{5}.*img,img])


%% Nicer display 
figure
subplot(2,3,1), imagesc(x_cs), colormap(gray(256)), title('Traditional CS')
subplot(2,3,2), imagesc(x_l1), colormap(gray(256)), title('L1-Analysis - db4')
subplot(2,3,3), imagesc(x_TV), colormap(gray(256)), title('TV Minimization')
subplot(2,3,4), imagesc(x_sos), colormap(gray(256)), title('db4 - Sum of squares')
subplot(2,3,5), imagesc(x_fcs), colormap(gray(256)), title('FCS - db4')
subplot(2,3,6), imagesc(x_fcs_TV), colormap(gray(256)), title('FCS - TV')



%% Save results and images
output_folder = './results/';
experiment_name = [experiment, '_N', num2str(n), '_wave', wave_name, '_nbSamples', sampling_ratio, '/'];

dir_path = [output_folder, experiment_name];

if isfolder(dir_path) == false
    mkdir(dir_path)
end

extensions = {'jpg', 'png', 'tif'};

for one_extension=1:length(extensions)
    imwrite((x_fcs-min(x_fcs(:)))/(max(x_fcs(:))-min(x_fcs(:))),[dir_path, 'FCS_results.', extensions{one_extension}])
    imwrite((x_fcs_TV-min(x_fcs_TV(:)))/(max(x_fcs_TV(:))-min(x_fcs_TV(:))),[dir_path, 'FCS_TV.', extensions{one_extension}])
    imwrite((x_l1-min(x_l1(:)))/(max(x_l1(:))-min(x_l1(:))),[dir_path, 'L1-Analysis.', extensions{one_extension}])
    imwrite((x_TV-min(x_TV(:)))/(max(x_TV(:))-min(x_TV(:))),[dir_path, 'TV-Min.', extensions{one_extension}])
    imwrite((x_sos-min(x_sos(:)))/(max(x_sos(:))-min(x_sos(:))),[dir_path, 'Sum-of-Squares.', extensions{one_extension}])
    imwrite((x_cs-min(x_cs(:)))/(max(x_cs(:))-min(x_cs(:))),[dir_path, 'Traditional_CS_results.', extensions{one_extension}])
    
    % Pointwise difference for CS
    normalized_x = (x_cs - min(min(x_cs)))/(max(max(x_cs)) - min(min(x_cs)));
    cur_err = img - normalized_x;
    cs_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'Traditional_CS_error.', extensions{one_extension}])
    % Pointwise difference for TV
    normalized_x = (x_TV-min(min(x_TV)))/(max(max(x_TV))-min(min(x_TV)));
    cur_err = img - normalized_x;
    tv_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'TV-min_error.', extensions{one_extension}])
    % Pointwise difference for l1 analysis
    normalized_x = (x_l1-min(min(x_l1)))/(max(max(x_l1))-min(min(x_l1)));
    cur_err = img - normalized_x;
    l1_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'L1-Analysis_error.', extensions{one_extension}])
    % Pointwise difference for fcs
    normalized_x = (x_fcs-min(min(x_fcs)))/(max(max(x_fcs))-min(min(x_fcs)));
    cur_err = img - normalized_x;
    fcs_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'FCS_error.', extensions{one_extension}])
    % Pointwise difference for fcs with TV
    normalized_x = (x_fcs_TV-min(min(x_fcs_TV)))/(max(max(x_fcs_TV))-min(min(x_fcs_TV)));
    cur_err = img - normalized_x;
    fcs_TV_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'FCS_TV_error.', extensions{one_extension}])
    % Pointwise difference for sos
    normalized_x = (x_sos-min(min(x_sos)))/(max(max(x_sos))-min(min(x_sos)));
    cur_err = img - normalized_x;
    sos_l2error = norm(cur_err)
    imwrite((cur_err - min(cur_err(:))) / (max(cur_err(:)) - min(cur_err(:)) ), [dir_path, 'Sum-of-Squares_error.', extensions{one_extension}])
end

snr_cs = psnr(x_cs,img)
snr_l1 = psnr(x_l1,img)
snr_tv = psnr(x_TV,img)
snr_sos = psnr(x_sos,img)
snr_fcs = psnr(x_fcs,img)
snr_fcs_TV = psnr(x_fcs_TV,img)
snr_fcs_TV2 = psnr(x_fcs_TV2,img)

%% Save .eps files
% Have to deal differently for .eps files

imagesc(x_fcs);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'FCS_results.eps'],'-deps2','-r300');  %# Print the figure


imagesc(x_cs);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'CS_results.eps'],'-deps2','-r300');  %# Print the figure



imagesc(x_fcs_TV);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'FCS_TV.eps'],'-deps2','-r300');  %# Print the figure


imagesc(x_l1);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'L1-Analysis.eps'],'-deps2','-r300');  %# Print the figure


imagesc(x_sos);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'Sum-of-Squares.eps'],'-deps2','-r300');  %# Print the figure


imagesc(x_TV);                     %# Plot the image
set(gca,'Units','normalized',...  %# Set some axes properties
        'Position',[0 0 1 1],...
        'Visible','off');
set(gcf,'Units','pixels',...      %# Set some figure properties
        'Position',[100 100 size(img,2) size(img,1)]);
print(gcf,[dir_path, 'TV-min.eps'],'-deps2','-r300');  %# Print the figure


%% Save parts of the workspace
save([dir_path, 'recovered_data.mat'], 'f_mask', 'x_fcs', 'x_TV', 'x_fcs_TV', 'x_cs', 'x_sos', 'x_l1')
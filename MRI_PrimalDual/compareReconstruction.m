function [denoised_outputs, times, img] = compareReconstruction(I, f_mask, masks, options)

if isfield(options, 'img_noise')
    img = I + options.img_noise*randn(size(I));
else
    img = I; % Makes sure the variable is defined when we need to use it. 
end

if isfield(options, 'noise_lvl')
    noise_lvl = options.noise_lvl; 
else
    noise_lvl = 0;
end

[m,n] = size(I);

%% Check do_XX
nb_doers = 9;
doers = {'do_cs', 'do_sos','do_fcs', 'do_ftv', 'do_l1synthesis', 'do_l1analysis', 'do_tv', 'do_fl1s', 'do_fl1a'};
for one_doer_idx=1:nb_doers
    cur_doer = doers{one_doer_idx};
    if ~isfield(options,cur_doer) % makes sure we haven't forgotten to set a certain field
        options = setfield(options, cur_doer, false);
    end
end

%% Generate the measurements

% The following should be part of the options
noise_type = options.noiseName;
% Remark: One should check that the do_XX have been set. Actually, we would
% default them to false, just to make sure, once we decide to implement the
% checks. (see box above)
%
% Single sensor processing? 
if options.do_cs || options.do_tv || options.do_l1analysis || options.do_l1synthesis
    switch noise_type
        case 'Normal' 
            img_f = c_fft_2d(img) + noise_lvl*randn(m,n);
        case 'Speckle' 
            img_f = imnoise(c_fft_2d(img), noise_type, noise_lvl); 
        case 'SnP'
            img_f = imnoise(c_fft_2d(img), 'salt & pepper', noise_lvl); 
    end
            
    y_rhs = f_mask.*img_f;
end
% 
% Fused sensor processing?
if options.do_fcs || options.do_ftv || options.do_fl1analysis || options.do_fl1synthesis || options.do_sos
    nb_channels = options.nb_channels;
    % Generate the multiple channel for the measurements.
    for one_channel=1:nb_channels
        
        switch noise_type
            case 'Normal'
                img_channels{one_channel} = c_fft_2d(masks{one_channel}.*img) + noise_lvl*randn(m,n);
            case 'Speckle'
                img_channels{one_channel} = imnoise(c_fft_2d(masks{one_channel}.*img), noise_type, noise_lvl);
            case 'SnP'
                img_channels{one_channel} = imnoise(c_fft_2d(masks{one_channel}.*img), 'salt & pepper', noise_lvl);
        end
        
        y_channel{one_channel} = f_mask.*img_channels{one_channel};
    end
    
end

%% Load the wavelet transforms
% Again, this should be check for existence of the fields if we were to do
% this nicely. 
WA_op = options.WA_op; % Wavelet analysis operator
WS_op = options.WS_op; % Wavelet synthesis operator

%% Let's do this! 
%% General Compressed sensing
if options.do_cs
    disp("Computing the recovery using traditional compressed sensing and DRS... ")
    
    % Call the primal dual algorithm
    disp('          ... Douglas-Rachford Splitting in progress ... ')
    prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    
    A=@(f) f_mask.*(c_fft_2d(f));
    At = @(f,gamma) real(f +c_ifft_2d(y_rhs-A(f)));
%     At = @(f,gamma) f +c_ifft_2d(y_rhs-A(f));
    [x, tout]= DouglasRachfordMin(zeros(size(y_rhs)),prox_l1_straight,At,options.DRoptions);
    denoised_outputs.cs = abs(x); times.cs = tout;
end

%% TV recovery 
if options.do_tv
    disp("Computing the recovery by total variation minimization... ")
    
    
    A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
    Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
        +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));
    
    abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
    sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
    prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
    prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
    
    prox_G = @(z) z + real(c_ifft_2d(y_rhs - f_mask.*c_fft_2d(z)));
    
    % Call the primal dual algorithm
    disp('          ... primal dual algorithm running ... ')
    [x, tout] = PrimalDualMin(zeros(m,n), prox_Fstar, prox_G, A, Astar, options.PDoptions);
    denoised_outputs.tv = x; times.tv = tout;
end

%% L1 synthesis 
if options.do_l1synthesis
    disp("Computing the recovery using l1 analysis... ")

    
    % Call the primal dual algorithm
    disp('          ... Douglas-Rachford Splitting algorithm running ... ')
    prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    
    A=@(f) f_mask.*(c_fft_2d(WS_op(f)));
    At = @(f,gamma) f +WA_op(real(c_ifft_2d(y_rhs-A(f))));
%     At = @(f,gamma) f +WA_op(c_ifft_2d(y_rhs-A(f)));
    
    [x, tout] = DouglasRachfordMin(zeros(size(At(y_rhs))),prox_l1_straight,At,options.DRoptions);
    denoised_outputs.l1s = abs(WS_op(x)); times.l1s = tout;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 analysis compressed sensing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.do_l1analysis
    disp("Computing the recovery using l1 analysis... ")
    
    % Call the DRS algorithm
    disp('          ... Douglas-Rachford Splitting algorithm running ... ')
    prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
    prox_l1_analysis = @(z,sigma) WS_op(prox_l1_straight(WA_op(z),sigma));
    
    A=@(f) f_mask.*(c_fft_2d(f));
    At = @(f,gamma) f +real(c_ifft_2d(y_rhs-A(f)));
%     At = @(f,gamma) f +c_ifft_2d(y_rhs-A(f));
    [x,tout]= DouglasRachfordMin(zeros(size(y_rhs)),prox_l1_analysis,At,options.DRoptions);
    denoised_outputs.l1a = abs(x); times.l1a = tout;
end


%% Fused compressed sensing
if options.do_fcs
    disp("Computing the recovery using fused CS sensing and DRS... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    fcs_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        % Call the primal dual algorithm
        disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
        A=@(f) f_mask.*(c_fft_2d(f));
        At = @(f,gamma) f +real(c_ifft_2d(y_rhs-A(f)));
%         At = @(f,gamma) f +c_ifft_2d(y_rhs-A(f));
        [x_local{one_channel}, tout] = DouglasRachfordMin(zeros(size(y_channel{one_channel})),prox_l1_straight,At,options.DRoptions);
        x = x + abs(x_local{one_channel});
        fcs_stops(one_channel) = tout;
    end
    
    tstart = tic;
    S = zeros(m,n);
    for i=1:nb_channels
        S = S + masks{i};
    end
    x = x./S;
    fcs_stops(nb_channels+1) = toc(tstart);
    
    denoised_outputs.fcs = x; times.fcs = fcs_stops;
    
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fused L1 analysis compressed sensing 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.do_fl1analysis
    disp("Computing the recovery using l1 analysis FCS ... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    fcs_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
        prox_l1_analysis = @(z,sigma) WS_op(prox_l1_straight(WA_op(z),sigma));
        
        A=@(f) f_mask.*(c_fft_2d(f));
        At = @(f,gamma) f +real(c_ifft_2d(y_channel{one_channel}-A(f)));
%         At = @(f,gamma) f +c_ifft_2d(y_channel{one_channel}-A(f));
        
        
        % Call the primal dual algorithm
        disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        [x_local{one_channel}, tout] = DouglasRachfordMin(zeros(size(At(y_channel{one_channel}))),prox_l1_analysis,At,options.DRoptions);
        x = x + abs(x_local{one_channel});
        fcs_stops(one_channel) = tout;
    end
    
    fcs_l1a_start = tic;
    S = zeros(m,n);
    for i=1:nb_channels
        S = S + masks{i};
    end
    x = x./S;
    fcs_stops(nb_channels+1) = toc(fcs_l1a_start);
    
    denoised_outputs.fl1a = x; times.l1a = fcs_stops;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fused L1 synthesis compressed sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.do_fl1synthesis
    disp("Computing the recovery using l1 synthesis FCS ... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    fcs_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
        
        A=@(f) f_mask.*(c_fft_2d(WS_op(f)));
        At = @(f,gamma) f +WA_op(real(c_ifft_2d(y_channel{one_channel}-A(f))));
%         At = @(f,gamma) f +WA_op(c_ifft_2d(y_channel{one_channel}-A(f)));
        
        
        % Call the primal dual algorithm
        disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        [x_local{one_channel}, tout] = DouglasRachfordMin(zeros(size(At(y_channel{one_channel}))),prox_l1_straight,At,options.DRoptions);
        x = x + abs(WS_op(x_local{one_channel}));
        fcs_stops(one_channel) = tout;
    end
    
    fcs_l1s_start = tic;
    S = zeros(m,n);
    for i=1:nb_channels
        S = S + masks{i};
    end
    x = x./S;
    fcs_stops(nb_channels+1) = toc(fcs_l1s_start);
    denoised_outputs.fl1s = x; times.fl1s = fcs_stops;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Faster Fused L1 synthesis compressed sensing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~options.do_fl1synthesis
    disp("Computing the recovery using l1 synthesis FCS ... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    fcs_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
        
        cur_mask = options.masks{one_channel};
        norm_of_col = sum(cur_mask.^2,1); % This is only for the sine masks. Basically, the projection operator should be adapted. 
        used_columns = (norm_of_col > .01);
        
        proj_small_set = @(x) x(:,used_columns);
        
        
        
        A=@(f) f_mask.*(c_fft_2d(back_to_big(WS_op(f), used_columns, m,n)  )  );
        At = @(f,gamma) f +WA_op(proj_small_set(c_ifft_2d(y_channel{one_channel}-A(f))));
        
%         size(At(y_channel{one_channel}))
        disp(['Hello there!'])
        size(proj_small_set(y_channel{one_channel}))
        disp(['Did we pass this guy?'])
        % zeros(size(proj_small_set(c_ifft_2d(y_channel{one_channel}))) )
        % Call the primal dual algorithm
        disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        [x_local{one_channel}, tout] = DouglasRachfordMin(zeros(size(proj_small_set(c_ifft_2d(y_channel{one_channel}))) ),prox_l1_straight,At,options.DRoptions);
        x = x + abs(WS_op(x_local{one_channel}));
        fcs_stops(one_channel) = tout;
    end
    
    fcs_l1s_start = tic;
    S = zeros(m,n);
    for i=1:nb_channels
        S = S + masks{i};
    end
    x = x./S;
    fcs_stops(nb_channels+1) = toc(fcs_l1s_start);
    denoised_outputs.fl1s = x; times.fl1s = fcs_stops;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fused TV mminimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.do_ftv
    disp("Computing the recovery by fused total variation minimization... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    fcs_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        % Proximal operator for the l2 data fidelity term
        A = @(f) [f - circshift(f,1,1) ,f - circshift(f,1,2)];
        Astar = @(f) (f(:, 1:n)- circshift(f(:, 1:n),-1,1) ) ...
            +(f(:,n+1:end)-circshift(f(:,n+1:end),-1,2));
        
        abs2 = @(z) repmat( sqrt(abs(z(:,1:n)).^2 +abs(z(:,1+n:end)).^2), 1,2);
        sign2 = @(z) [real(sign(z(:,1:n)+1i*z(:,1+n:end))), imag(sign(z(:,1:n)+1i*z(:,1+n:end)))];
        prox_F = @(sigma, z) sign2(z).*(abs2(z)-sigma).*(abs2(z)>=sigma) ;
        prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);
        
        prox_G = @(z) z + real(c_ifft_2d(y_channel{one_channel}- f_mask.*c_fft_2d(z)));
        
        
        % Call the primal dual algorithm
        disp(['          ... primal dual algorithm running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        [x_local{one_channel}, tout] = PrimalDualMin(zeros(m,n), prox_Fstar, prox_G, A, Astar, options.PDoptions);
        x = x + x_local{one_channel};
        fcs_stops(one_channel) = tout;
    end
    
    fcs_TV_start = tic;
    S = zeros(m,n);
    for i=1:nb_channels
        S = S + masks{i};
    end
    x = x./S;
    fcs_stops(nb_channels+1) = toc(fcs_TV_start);
    denoised_outputs.ftv = x; times.ftv = fcs_stops;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L1 analysis compressed sensing -- sum of squares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if options.do_sos
    disp("Computing the recovery using l1 analysis SoS ... ")
    
    % Initialize the end result with zeros
    x = zeros(m,n);
    
    sos_stops = zeros(nb_channels+1,1);
    for one_channel=1:nb_channels
        prox_l1_straight = @(z, sigma) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
        prox_l1_analysis = @(z,sigma) WS_op(prox_l1_straight(WA_op(z),sigma));
        
        A=@(f) f_mask.*(c_fft_2d(f));
        At = @(f,gamma) f +real(c_ifft_2d(y_channel{one_channel}-A(f)));
%         At = @(f,gamma) f +c_ifft_2d(y_channel{one_channel}-A(f));
        
        
        % Call the primal dual algorithm
        disp(['          ... DRS running ... Channel number ', num2str(one_channel), ' from ', num2str(nb_channels)])
        [x_local{one_channel}, tout] = DouglasRachfordMin(zeros(size(y_channel{one_channel})),prox_l1_analysis,At,options.DRoptions);
        x = x + abs(x_local{one_channel}).^2;
        sos_stops(one_channel) = tout;
    end

    sos_start = tic;
    x = sqrt(x);
    sos_stops(nb_channels+1) = toc(sos_start);
    
    denoised_outputs.sos = x; times.sos = sos_stops;
end



%% Auxiliary functions
function backProjection = back_to_big(cur_data, cur_columns, m,n)
backProjection = zeros(m,n);
backProjection(:,cur_columns) = cur_data;
figure
imagesc(cur_data), colormap(gray(256))
figure
imagesc(backProjection), colormap(gray(256))
return
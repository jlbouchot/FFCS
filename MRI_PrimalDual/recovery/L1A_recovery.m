function x_hat = L1A_recovery(y_rhs,f_mask, x0, measurement, measurementStar, PDoptions, WAVEoptions)

if nargin < 4 || isempty(measurement)
    measurement = @(x) c_fft_2d(x);
end

if nargin < 5 || isempty(measurementStar)
    measurementStar = @(f) c_ifft_2d(f);
end

if nargin < 3 || isempty(x0)
    x0 = zeros(size(measurementStar(y_rhs)));
end
disp('L1Analysis recovery')
wave_name = WAVEoptions.wname;
nb_level = WAVEoptions.nblvl;

% Get the structure first, this will not change through time
[~,s]=wavedec2(x0,nb_level,wave_name);

% Analysis operator
A = @(f) wavedec2(f,nb_level,wave_name);
% Synthesis operator
Astar = @(f) waverec2(f,s,wave_name);

prox_F = @(sigma, z) sign(z).*(abs(z)-sigma).*(abs(z)>=sigma) ;
prox_Fstar = @(sigma,z) z-sigma*prox_F(1/sigma, z/sigma);

prox_G = @(z) z + real(measurementStar(y_rhs- f_mask.*measurement(z)));

% Call the primal dual algorithm
disp('          ... primal dual algorithm running ... ')
x_hat = PrimalDualMin(x0, prox_Fstar, prox_G, A, Astar, PDoptions);
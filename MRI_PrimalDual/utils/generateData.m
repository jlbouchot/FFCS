function [img, experiment, m,n] = generateData(which_data, m)
if nargin < 1 || isempty(which_data)
    which_data = 1; 
end

if nargin < 2 || isempty(m)
    m = 256;
end

if which_data == 1
    n = m;
    img = phantom(n); % Shep-Logan CT Phantom
    
    img = img/max(img(:));
    
    experiment = 'artificial'; % used for automatic saving
elseif which_data == 2
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
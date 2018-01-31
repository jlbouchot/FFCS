function xhat = l1minl2data(A, y, noise,D)
% L1MINL2DATA Computes the l1 analysis recovery with CVX
% 
%   * A is the m-by-N sensing matrix
%   * y is the m-dimensional measurement vector (y = Ax + e)
%   * noise is a l2-bound on the noise ||e||_2 <= noise
%   * D corresponds to a N-by-p dictionary such that D^Tx is sparse
% 
% 
%   Author  : Jean-Luc Bouchot
%   Contact : bouchot@mathc.rwth-aachen.de
%   This is for scholarly uses only and should be acknowledged if used.
% 
%   Created on January 31st 2018. 
%   Last modified on January 31st 2018. 

% Obviously we can still improve here and there! 

N = size(A,2);

cvx_begin
    variable z(N)
    minimize norm(D'*z,1)
    norm(A*z-y,2) <= noise;
cvx_end

xhat = z;
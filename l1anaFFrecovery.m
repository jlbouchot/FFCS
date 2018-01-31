function xhat = l1anaFFrecovery(A, measures, noises, Pi, Di, Sinv)
% L1ANAFFRECOVERY Computes the l1 analysis recovery with fusion frames.
%
%   XHAT = L1ANAFFRECOVERY(A, measures, noises, Pi, Di, S) returns the
%   vector obtained via compressed sensing-fusion frame recovery. 
%   
%   * DI is a cell structure containing the n sparsifying dictionaries (n
%   being the number of projections in our fusion frame structure). In case
%   a unique dictionnary sparsifies all the projections, Di = D is a unique
%   matrix. 
%   * MEASURES is an m-by-n matrix containing (column wise) the n
%   m-dimensional measurement vectors taken (separately) at the n stations
%   * PI corresponds to the projection matrices. These are represented by
%   n cells each containing a matrix. 
%   A is the (unique for nw) m-vy-N sensing matrix used at every station.
%   Hence the measures are expected to be measures(:,j) = A*Pi{j} plus some
%   noise. 
%   NOISES is an n-dimensional vector of  either 
%       * floats: in this case we will use classical l1 analysis
%       minimization
%       * integers: this corresponds to the case where we would use an
%       iterative algorithm such as HTP
%   SINV corresponds to the inverse (fusion) frame operator (Note: in a
%   near future, it would be good to have some a way to use the frame
%   operator in case Sinv is not known)
%   
%   TODOs in the future:
%       * Have a switching variable for the choice of local recovery algo
%       * Have a smarter way to deal with the inverse frame operator (Sinv
%       or the frame algorithm?)
%       * Tolerate different A for each measurement vector: A{j}*Pi{j}.
%
%   Authors : Roza Aceska, Jean-Luc Bouchot, Shidong Li
%   Contact : bouchot@mathc.rwth-aachen.de
%   This is for scholarly uses only and the following papers should be
%   acknowledged / cited in case of use: 
%   * Roza Aceska, Jean-Luc Bouchot, and Shidong Li, "Local sparsity and
%   recovery of fusion frame structured signals", Circuits, Systems &
%   Signal Processing, 2018.
%   * Roza Aceska, Jean-Luc Bouchot, and Shidong Li, "Fusion Frames and
%   Distributed Sparsity", Contemporary Mathematics, 2018. 

% Should add some safeguards to deal with stupid users not entering the
% correct arguments.

nbStations = size(measures,2); 
N = size(A,2);

xhat = zeros(N,1);

for one_station=1:nbStations % In a neat fast implementation, this is done in parallel
    % Basically solve the following problem:
    % minimize ||Di{one_station}z||_1, 
    % s.t. ||APi{one_station}z - y(:,one_station)||_2 \leq noises(station)
    z = l1minl2data(A*Pi{one_station}, measures(:,one_station), noises(one_station), Di{one_station}); % (A, y, noise,D)
%     z
%     cvx_begin
%         variable z(N)
%         minimize norm(Di{one_station}'*z,1)
%         norm(A*Pi{one_station}*z-measures(:,one_station),2) <= noises(one_station);
%     cvx_end
%     z
    xhat = xhat + z;
end


xhat = Sinv*xhat;

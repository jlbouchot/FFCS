# Recovery of Fusion Frame Structured Signal via Compressed Sensing 

This repository contains some accompanying files used for the paper 

Roza Aceska, Jean-Luc Bouchot, and Shidong Li,
**Local sparsity and recovery of fusion frames structured signals**
_arXiv:1604.00424_
and it should be cited anywhere these files are used together with 
Roza Aceska, Jean-Luc Bouchot, and Shidong Li,
**Fusion Frames and Distributed Sparsity**
__Contemporary Mathematics__, 2018. 

Note that the example requires the use of a full Haar matrix. This was just for sake of simplicity of the implementation. We used the file obtained [here](http://de.mathworks.com/matlabcentral/fileexchange/4619-haarmtx?focused=5054050&tab=function "To the Haar matrix!"). Any implementation of your liking should work just fine, as soon as it has the name `haarmtx`.

The l1-analysis optimization is done using [CVX](http://cvxr.com/cvx/download/ "Convex optimization without sweating").

We welcome any questions / remarks / complaints / improvements sent to 
bouchot@mathc.rwth-aachen.de .
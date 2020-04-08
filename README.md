# Recovery of Fusion Frame Structured Signal via Compressed Sensing 

This repository contains some accompanying files used for the papers


Roza Aceska, Jean-Luc Bouchot, and Shidong Li,
**Local sparsity and recovery of fusion frames structured signals**,
_arXiv:1604.00424_ , (or [here](./Papers/ABL_SigProc.pdf) for the pdf file)
and it should be cited anywhere these files are used together with 
Roza Aceska, Jean-Luc Bouchot, and Shidong Li,
**Fusion Frames and Distributed Sparsity**,
__Contemporary Mathematics__, 2018, which is also available [here as pdf](./Papers/ABL_conm.pdf).

# Running the experiments

They are two main experiments in this repository which should be working independently from one another. 
All files are implemented using Matlab.

## Denoising of Doppler signals. 

The first one deals with denoising of Doppler signals
This example requires the use of a full Haar matrix. One could easily implement a version without the need of actually implementing it, which is especially bad once the dimensionality increases, but it is not the point of this reposotory and would only hinder the understanding of the proposed approach. 
We used the file obtained [here](http://de.mathworks.com/matlabcentral/fileexchange/4619-haarmtx?focused=5054050&tab=function "To the Haar matrix!"). 
Any implementation of your liking should work just fine, as soon as it has the name `haarmtx`. 
Or you can simply fiddle with the names.

The l1-analysis optimization is done using [CVX](http://cvxr.com/cvx/download/ "Convex optimization without sweating").

## MR images reconstruction

The bigger part of the demonstration deals with denoising of medical images, specifically in the context of Magnetic Resonnance Imaing. 

### Overview

It is known that such images enjoy a sparse representation in wavelets and / or a minimal total variation. 
Both approaches are implemented in the [MRI_PrimalDual](MRI_PrimalDual repository). 
This folder contains all the required files (except some signal processing files, see  the next section) among which you may find (but not limited to)
* A Primal dual algorithm, freshly implemented. It is implemented in a way that should handle very high dimensional data. It has been tested with a 2048*2048 MRI (i.e. over 4 Millions of parameters)
* A Douglas-Rachford optimizer. 
Both optimization functions are using proximal methods for the recovery of sparse signals. 
They are implemented assuming equality contraints in the data fidelity terms but could, potentially, be adapted to the inequality case. 

To play around, simply download the folder and run the `LastTestsForRevision.m`. 

### Some links
* (Link to the numerical-tours)[https://www.numerical-tours.com/matlab/], needed to run the script. It contains a bunch of functions for signal processing and analysis
* (Link to Chun and Adcock's paper)[https://ieeexplore.ieee.org/document/7917266/] dealing with parallel acquisition
* (Link to the variable density sampling schemes)[https://arxiv.org/abs/1311.6039] by Chauffert and co-authors. 

## Contacts / Remarks
We welcome any questions / remarks / complaints / improvements sent to 
jlbouchot@gmail.com .

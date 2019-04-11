function [Peaks, RefinedPeaks]=ColumnFinderSeries(Series, ProjPeaksGauss, GaussWindow, Guess, lb, ub, Noise, RoseCriterion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %%%                                 INPUTS
%
% Series - A time series of images of a set of particles or atoms
% RefinedLattice - A list of lattice points around which you want to look
%                  for atom positions. 
% PeakfindWindow - The diameter of circle around each lattice point over which you want
%                  to look for local maxima
% CentroidWindow - The diameter of circle over which you want to do a centroid fit
%                  about the local maximum. 
% GaussWindow    - The width of a box around the centroid position within
%                  which you will fit a 2D gaussian. 
% Noiselevel     - A threshold below which a peak will not be fitted.
%
% %%%                              DESCRIPTION 
%
% peakrefiner takes an image series (Series) with dimensions nx,ny,nz.
% It first runs locmax on each image, to find maxima within a given radius
% of a lattice point
% It then runs cntrd, and then fitgauss2d on each image to refine the position of the atomic
% columns with sub-pixel accuracy by fitting a centroid and then a 2D gaussian to the data 
 
% %%%                                 OUTPUTs
%
% The first output PeaksGauss, is a 3 dimensional array. with dimensions (number of peaks detected, 5, nz)
% Each 2 - dimensional layer of Peaks contains 6 columns of data related to each image in the series A. 
% In order from 1 to 6, these are:
% 1) the amplitude of each gaussian fitted to the data, 
% 2) the x co-ordinates of the center each gaussian, 
% 3) the x-variance of each gaussian
% 4) the y co-ordinates of the center each gaussian, 
% 5) the y-variance of each gaussian
% 6) An I.D. number associated with each atom/particle.
%
% The second output, RefinedPeaksGauss is a 2 dimensional array with 7 columns.
% This is essentially just a reshaped version of peaks. The first 6 columns
% contain the same information as in Peaks, but for all frames. The 7th
% column lists the frame number for each data point. This second format is
% useful for making certain plots (like the tracks snd rmsds). 
%
% The scripts locmax, and cntrd are from the MATLAB particle tracking package by
% Daniel Blair and Eric Dufresne, which is based on the IDL particle tracking software
% by David Grier, John Crocker, and Eric Weeks. 
%
% References J. C. Crocker, D. G. Grier, Methods of Digital Video 
% Microscopy for Colloidal Studies. J. Colloid Interface Sci.179, 298–310 
% (1996). doi:10.1006/jcis.1996.0217.
% 
% D. Blair, E.Dufresne, (The Matlab Particle Tracking Code Repository, at 
% http://physics.georgetown.edu/matlab/)
%
% Gaussian fitting is performed using a modified version of the
% D2GaussFunction code written by Gero Nootz https://www.mathworks.com/matlabcentral/fileexchange/37087-fit-2d-gaussian-function-to-data
%
% By Barnaby Levin, ASU, 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Asize=size(Series);
sizeAsize = size(Asize);

if sizeAsize(2)==2
%     MFSeries=medfilt2(Series,[m,n]);
    Gauss = iterategauss2d2(Series,ProjPeaksGauss,GaussWindow, Guess, lb, ub, Noise, RoseCriterion);
    Peaks(:,:) = Gauss(:,:);
    sizeP = size(Peaks);
    RefinedPeaks = zeros(sizeP(1), sizeP(2)+1);
    RefinedPeaks(:, 1:sizeP(2))=Peaks(:,:);
    RefinedPeaks(:, sizeP(2)+1)=1;
    
else
    numim = Asize(3);
    Data = zeros(1000,11,numim); 
    
    for n = 1:numim
       Gauss = iterategauss2d2(Series(:,:,n),ProjPeaksGauss,GaussWindow, Guess, lb, ub, Noise, RoseCriterion);
       Datasize = size(Gauss);
       Data(1:Datasize(1),:,n)= Gauss(:,:);

    end

    Peaks = Data(1:Datasize, :, :);
    sizeP = size(Peaks);
    RefinedPeaks = zeros(sizeP(1)*numim, sizeP(2)+1);
    for m = 1:numim
        RefinedPeaks(1+sizeP(1)*(m-1):sizeP(1)*m, 1:sizeP(2))=Peaks(:,:,m);
        RefinedPeaks(1+sizeP(1)*(m-1):sizeP(1)*m, sizeP(2)+1)=m;
    end
        
end
end
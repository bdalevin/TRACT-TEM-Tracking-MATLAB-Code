function Jumps=track_gaussjumps2(TrackedPeaks, scale)
% This function takes as its input a 3 dimensional array (the output of 
% track_descrambler.m), and calculates the root mean squared displacement 
% of atoms.
% The second input is the size of 1 pixel in the image in nm
% The output is a 2D array, listing the mean x and y coordinates of atomic 
% columns, along with 

sizeJ = size(TrackedJumpsGauss);
% Pre-allocate memory
Jumps = zeros (sizeJ(3),9);
coords = zeros(sizeJ(1),2);
for n = 1:sizeJ(3)

    % Work out mean coordinates for origin
    coords=TrackedJumpsGauss(:,2:2:4,n);
    %coords( ~any(coords,2), : ) = []; % Remove any zeros 
    meanx = mean(coords(:,1));
    meany = mean(coords(:,2));
    sizecoords = size(coords);
    
    % Work out coordinates relative to mean
    coords0 = zeros(sizecoords(1),5);
    coords0(:,1) = coords(:,1)-meanx;
    coords0(:,2) = coords(:,2)-meany;
    coords0(:,3) = TrackedJumpsGauss(:,7,n);
    coords0(:,4) = TrackedJumpsGauss(:,8,n);
    
    
    % Work out root squared displacement for every atom
    coords0(:,5) = sqrt(coords0(:,3).^2+coords0(:,4).^2)*scale;
    
    % Work out Number of Jumps
    J = sum(coords0(:,5)>10);
    
    
    % Assign values to F
    Jumps(n,1) = meanx;
    Jumps(n,2) = meany;
    Jumps(n,3) = J;
    Jumps(n,4) = -0.025*log(J/(10^13)); % Activation Energy.
    
end

end


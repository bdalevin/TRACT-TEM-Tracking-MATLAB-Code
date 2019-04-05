function RMSDs=track_gaussrmsd2(TrackedPeaks, ProjPeaksGauss, scale)
% This function takes as its input a 3 dimensional array (the output of 
% track_descrambler.m), and calculates the root mean squared displacement 
% of atoms.
% The second input is the size of 1 pixel in the image in nm
% The output is a 2D array, listing the mean x and y coordinates of atomic 
% columns, along with 

sizeP = size(TrackedPeaks);
sizePPG=size(ProjPeaksGauss);
PPG = zeros(1,sizePPG(2),sizePPG(1));
for m = 1:sizePPG(1)
   PPG(1,:,m)= ProjPeaksGauss(m,:);
end

% Pre-allocate memory
RMSDs = zeros (sizeP(3),10);
coords = zeros(sizeP(1),2);
for n = 1:sizeP(3)

    % Work out mean coordinates for origin
    coords=TrackedPeaks(:,2:2:4,n);
    coords( ~any(coords,2), : ) = []; % Remove any zeros 
    meanx = PPG(1,2,n);
    meany = PPG(1,4,n);
    sizecoords = size(coords);
    
    % Work out coordinates relative to mean
    coords0 = zeros(sizecoords(1),3);
    coords0(:,1) = coords(:,1)-meanx;
    coords0(:,2) = coords(:,2)-meany;
    
    % Work out root squared displacement for every atom
    coords0(:,3) = sqrt(coords0(:,1).^2+coords0(:,2).^2);
    
    % Work out rmsd
    rmsd = mean(coords0(:,3));
    mrsd = median(coords0(:,3));
    rmsx = mean(sqrt(coords0(:,1).^2));
    rmsy = mean(sqrt(coords0(:,2).^2));
    
    % Assign values to F
    RMSDs(n,1) = meanx;
    RMSDs(n,2) = meany;
    RMSDs(n,3) = rmsd;
    RMSDs(n,4) = rmsd*scale*1000; % Rmsd in picometers.
    RMSDs(n,5) = rmsx;
    RMSDs(n,6) = rmsx*scale*1000; % Rmsx in picometers.
    RMSDs(n,7) = rmsy;
    RMSDs(n,8) = rmsy*scale*1000; % Rmsy in picometers.
    RMSDs(n,9) = TrackedPeaks(1,9,n);
    RMSDs(n,10) = mrsd;
    
    
end

end


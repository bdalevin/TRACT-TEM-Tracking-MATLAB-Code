function [TrackedJumps, Jumps] = track_jumps(TrackedPeaks, scale, Threshold)

SizeTP = size(TrackedPeaks);
TrackedJumps=zeros(SizeTP(1),SizeTP(2)+3,SizeTP(3)); 
TrackedJumps(:,1:SizeTP(2),:)=TrackedPeaks(:,:,:);

TrackedJumps(2:SizeTP(1),SizeTP(2)+1,:)= TrackedPeaks(2:SizeTP(1),1,:)-TrackedPeaks(1:SizeTP(1)-1,1,:);
TrackedJumps(2:SizeTP(1),SizeTP(2)+2,:)= TrackedPeaks(2:SizeTP(1),2,:)-TrackedPeaks(1:SizeTP(1)-1,2,:);
TrackedJumps(:,SizeTP(2)+3,:)=sqrt((TrackedJumps(:,SizeTP(2)+1,:).*TrackedJumps(:,SizeTP(2)+1,:))+((TrackedJumps(:,SizeTP(2)+2,:).*TrackedJumps(:,SizeTP(2)+2,:))));

sizeJ = size(TrackedJumps);
% Pre-allocate memory
Jumps = zeros (sizeJ(3),4);
coords = zeros(sizeJ(1),2);
for n = 1:sizeJ(3)

    % Work out mean coordinates for origin
    coords=TrackedJumps(:,1:2,n);
    %coords( ~any(coords,2), : ) = []; % Remove any zeros 
    meanx = mean(coords(:,1));
    meany = mean(coords(:,2));
    sizecoords = size(coords);
    
    % Work out coordinates relative to mean
    coords0 = zeros(sizecoords(1),5);
    coords0(:,1) = coords(:,1)-meanx;
    coords0(:,2) = coords(:,2)-meany;
    coords0(:,3) = TrackedJumps(:,SizeTP(2)+1,n);
    coords0(:,4) = TrackedJumps(:,SizeTP(2)+2,n);
    
    
    % Work out root squared displacement for every atom
    coords0(:,5) = sqrt(coords0(:,3).^2+coords0(:,4).^2)*1000*scale;
    
    % Work out Number of Jumps
    J = sum(coords0(:,5)>Threshold);
    
    % Assign values to F
    Jumps(n,1) = meanx;
    Jumps(n,2) = meany;
    Jumps(n,3) = J; % Number of times jump over threshold
    if J ==0
        Jumps(n,4)=0.75;
    else
    Jumps(n,4) = -0.025*log(J/(10^13)); 
    end
end
end
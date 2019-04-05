function [TrackedJumpsGauss, Jumps] = track_gaussjumps(TrackedPeaksGauss, scale, Threshold)

SizeTPG = size(TrackedPeaksGauss);
TrackedJumpsGauss=zeros(SizeTPG(1),SizeTPG(2)+3,SizeTPG(3)); 
TrackedJumpsGauss(:,1:SizeTPG(2),:)=TrackedPeaksGauss(:,:,:);

TrackedJumpsGauss(2:SizeTPG(1),SizeTPG(2)+1,:)= TrackedPeaksGauss(2:SizeTPG(1),2,:)-TrackedPeaksGauss(1:SizeTPG(1)-1,2,:);
TrackedJumpsGauss(2:SizeTPG(1),SizeTPG(2)+2,:)= TrackedPeaksGauss(2:SizeTPG(1),4,:)-TrackedPeaksGauss(1:SizeTPG(1)-1,4,:);
TrackedJumpsGauss(:,SizeTPG(2)+3,:)=sqrt((TrackedJumpsGauss(:,SizeTPG(2)+1,:).*TrackedJumpsGauss(:,SizeTPG(2)+1,:))+((TrackedJumpsGauss(:,SizeTPG(2)+2,:).*TrackedJumpsGauss(:,SizeTPG(2)+2,:))));

sizeJ = size(TrackedJumpsGauss);
% Pre-allocate memory
Jumps = zeros (sizeJ(3),5);
coords = zeros(sizeJ(1),2);
Freq = 10^13;
for n = 1:sizeJ(3)

    % Work out mean coordinates for origin
    coords=TrackedJumpsGauss(:,2:2:4,n);
    
    for k = 1:round(SizeTPG(1)/2)
        if isnan(coords(k,1))
            Index=find(~isnan(coords(k:SizeTPG(1),1)));
            coords(k,:)=coords(Index(1)+k-1,:);
        end 
    end
    for m = round(SizeTPG(1)/2):SizeTPG(1)
        if isnan(coords(m,1))
            Index=find(~isnan(coords(1:m,1)));
            Index=flipud(Index);
            coords(m,:)=coords(Index(1),:);
        end 
    end
    
    sizecoords = size(coords);
    meanx = mean(coords(:,1));
    meany = mean(coords(:,2));
    
    
    % Work out coordinates relative to mean
    coords0 = zeros(sizecoords(1),5);
    coords0(:,1) = coords(:,1)-meanx;
    coords0(:,2) = coords(:,2)-meany;
    coords0(:,3) = TrackedJumpsGauss(:,SizeTPG(2)+1,n);
    coords0(:,4) = TrackedJumpsGauss(:,SizeTPG(2)+2,n);
    
    
    % Work out root squared displacement for every atom
    coords0(:,5) = sqrt(coords0(:,3).^2+coords0(:,4).^2)*1000*scale; %
    coords0(:,6) = sqrt(coords0(:,1).^2+coords0(:,2).^2)*1000*scale;
    
%     % Work out Number of Jumps
%     J = sum(coords0(:,6)>Threshold); % Number of jumps away from mean greater than threshold
    
    Counter = zeros(sizeJ(1),1);
    for num_frames = 2:sizeJ(1)
        if abs(coords0(num_frames,5)-coords0(num_frames-1,5))>Threshold && coords0(num_frames,6)>Threshold
            Counter(num_frames,1) = 1;
        else 
            Counter(num_frames,1) = 0;
        end
    end
    J = sum(Counter(:,1)); % Number of jumps away from previous frame greater than threshold
    
    
    
    % Assign values to F
    Jumps(n,1) = meanx;
    Jumps(n,2) = meany;
    Jumps(n,3) = J; % Number of times jump over threshold
    if J ==0
        Jumps(n,4)=-0.025*log(1/(Freq)); % Set to max energy if no motion
    else
    Jumps(n,4) = -0.025*log(J/(Freq)); 
    end
    Jumps(n,5) = n;
end
end
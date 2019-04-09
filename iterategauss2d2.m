function out=iterategauss2d2(im,peaks,sz, Guess, lb, ub, Noise, RoseCriterion)
% out=cntrd(im,mx,sz,interactive)
% 
% PURPOSE:  calculates the centroid of bright spots to sub-pixel accuracy.
%  Inspired by Grier & Crocker's feature for IDL, but greatly simplified and optimized
%  for matlab
% 
% INPUT:
% im: image to process, particle should be bright spots on dark background with little noise
%   ofen an bandpass filtered brightfield image or a nice fluorescent image
%
% mx: locations of local maxima to pixel-level accuracy from pkfnd.m
%
% sz: diamter of the window over which to average to calculate the centroid.  
%     should be big enough
%     to capture the whole particle but not so big that it captures others.  
%     if initial guess of center (from pkfnd) is far from the centroid, the
%     window will need to be larger than the particle size.  RECOMMENDED
%     size is the long lengthscale used in bpass plus 2.
%     
%
% interactive:  OPTIONAL INPUT set this variable to one and it will show you the image used to calculate  
%    each centroid, the pixel-level peak and the centroid
%
% NOTE:
%  - if pkfnd, and cntrd return more then one location per particle then
%  you should try to filter your input more carefully.  If you still get
%  more than one peak for particle, use the optional sz parameter in pkfnd
%  - If you want sub-pixel accuracy, you need to have a lot of pixels in your window (sz>>1). 
%    To check for pixel bias, plot a histogram of the fractional parts of the resulting locations
%  - It is HIGHLY recommended to run in interactive mode to adjust the parameters before you
%    analyze a bunch of images.
%
% OUTPUT:  a N x 4 array containing, x, y and brightness for each feature
%           out(:,1) is the x-coordinates
%           out(:,2) is the y-coordinates
%           out(:,3) is the brightnesses
%           out(:,4) is the square of the radius of gyration
%
% CREATED: Eric R. Dufresne, Yale University, Feb 4 2005
%  5/2005 inputs diamter instead of radius
%  Modifications:
%  D.B. (6/05) Added code from imdist/dist to make this stand alone.
%  ERD (6/05) Increased frame of reject locations around edge to 1.5*sz
%  ERD 6/2005  By popular demand, 1. altered input to be formatted in x,y
%  space instead of row, column space  2. added forth column of output,
%  rg^2
%  ERD 8/05  Outputs had been shifted by [0.5,0.5] pixels.  No more!
%  ERD 8/24/05  Woops!  That last one was a red herring.  The real problem
%  is the "ringing" from the output of bpass.  I fixed bpass (see note),
%  and no longer need this kludge.  Also, made it quite nice if mx=[];
%  ERD 6/06  Added size and brightness output ot interactive mode.  Also 
%   fixed bug in calculation of rg^2
%  JWM 6/07  Small corrections to documentation 

% if nargin==3
%    interactive=0;
%    Constraint = 10000000;
% end
% 
% if nargin==4
%    interactive=0;
% end

if sz/2 ~= floor(sz/2)
sz=sz+1;
end

if isempty(peaks)
    warning('there were no positions inputted into cntrd. check your pkfnd theshold')
    out=[];
    return;
end


r = sz/2; % Defines radius of window
% Separate IDs from coordinates
ID = peaks(:,9);
mx = peaks(:,2:2:4);
sig = peaks(:,3:2:5);
[nr,nc]=size(im);
%remove all potential locations within distance sz from edges of image
ind=find(mx(:,2) > r & mx(:,2) < nr-r);
mx=mx(ind,:);
ind=find(mx(:,1) > r & mx(:,1) < nc-r);
mx=mx(ind,:);

[nmx,crap] = size(mx);

%inside of the window, assign an x and y coordinate for each pixel
xl=zeros(2*r,2*r); % Defines window (diameter 2r)
for i=1:2*r
    xl(i,:)=(1:2*r); 
end
yl=xl';
GaussParams = zeros(nmx,9);
pts= zeros(nmx,9);
mx2 = round(mx);
%loop through all of the candidate positions
for i=1:nmx
    %Crop out a box around each maximum. Fit gaussian inside the box. Find position of centre of Gaussian 
    tmp=im((mx2(i,2)-r:mx2(i,2)+r),(mx2(i,1)-r:mx2(i,1)+r));
%     x0 = r; y0 = r; 
%     sigmax = 2; sigmay = 2;
%     Amp = max(max(tmp))-median(min(tmp));
%     GaussParams(i,:) = fit2dgauss(tmp,Amp,x0,sigmax,y0,sigmay, lb, ub, Scale, Noise);
    GaussParams(i,:) = fit2dgauss(tmp, Guess, lb, ub, Noise, RoseCriterion);

%     xavg = GaussParams(i,2)+r;
%     yavg = GaussParams(i,4)+r;

%     tmp1=im((mx(i,2)-round(0.6*r):mx(i,2)+round(0.6*r)),(mx(i,1)-round(0.6*r):mx(i,1)+round(0.6*r)));
%     tmp2=im((mx(i,2)-round(1.6*r):mx(i,2)+round(1.6*r)),(mx(i,1)-round(1.6*r):mx(i,1)+round(1.6*r)));
%     szt1=size(tmp1);
%     szt2=size(tmp2);
%     tmp1r=reshape(tmp1,szt1(1)*szt1(2),1);
%     tmp2r=reshape(tmp2,szt2(1)*szt2(2),1);
%     Med1=median(tmp1r);
%     Med2=median(tmp2r);
%     
%     if Med1<Constraint*Med2
%         GaussParams(i,2)=NaN;
%         GaussParams(i,4)=NaN;
%     end
    
%     if GaussParams(i,2)>Constraint*sig(i,1)
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     elseif  GaussParams(i,4)>Constraint*sig(i,2)
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     end

%     if GaussParams(i,3)>sig(i,1)
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     elseif  GaussParams(i,5)>sig(i,2)
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     end
    
%     % Set a threshold below which peak will not be indexed.
%     sztmp = size(tmp);
%     tmparray = reshape(tmp,sztmp(1)*sztmp(2),1);
%     y = std(tmparray);
%     if y < noiselevel
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     end
%     if GaussParams(i,1) < noiselevel
%             GaussParams(i,2)=0;
%             GaussParams(i,4)=0;
%     end
    
    % Convert co-ordinates from those of the box to those of the whole image, and save. 
    GaussParams(i,2) = mx(i,1)+(GaussParams(i,2))+0.5; % You need to add the plus 0.5 here because MATLAB starts counting arrays from 1. 
    GaussParams(i,4) = mx(i,2)+(GaussParams(i,4))+0.5; % Without the plus 0.5, everything gets plotted one pixel up, and one pixel left of where it should be.   
    pts(i,1:8)= GaussParams(i,1:8);
    pts(i,9) = ID(i);
    
    
end
out=pts;


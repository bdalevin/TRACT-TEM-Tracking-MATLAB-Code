function GaussParams = fit2dgauss(Image, Guess, lb, ub, Noise, RoseCriterion)
%fit2dgauss(Image,Amp,x1,sigmax,y1,sigmay, lb, ub, Scale, Noise)
% Fit a 2D gaussian function to data
%% PURPOSE:  Fit a 2D gaussian centroid to simulated data
% Uses lsqcurvefit to fit
%
% INPUT:
% 
%   MdataSize: Size of nxn data matrix
%   x0 = [Amp,x0,wx,y0,wy,fi]: Inital guess parameters
%   x = [Amp,x0,wx,y0,wy,fi]: simulated centroid parameters
%   noise: noise in % of centroid peak value (x(1)
%   InterpolationMethod = 'nearest' or 'linear' or 'spline' or 'cubic'
%       used to calculate cross-section along minor and major axis
%     
%
%
% NOTE:
%   The initial values in x0 must be close to x in order for the fit
%   to converge to the values of x (especially if noise is added)
%
% OUTPUT:  non
%
% CREATED: G. Nootz  May 2012
% 
%  Modifications:
%  non
%% ---------User Input---------------------

Imsize = size(Image); % Size of image
MdataSize = Imsize(1);
x = [-MdataSize/2:MdataSize/2-1]; y = [-MdataSize/2:MdataSize/2-1];
% I = Image;
% I(MdataSize/5:4*MdataSize/5,MdataSize/5:4*MdataSize/5);
% Z = Image-median(median(I));
% XAx(:,1) = Image(round(MdataSize/2),:)';
% YAx(:,1) = Image(:,round(MdataSize/2))';
% Diag1(:,1)=zeros(MdataSize,1);
% Diag2(:,1)=zeros(MdataSize,1);
% for pixval=1:MdataSize
%     Diag1(pixval,1)=Image(pixval,pixval);
%     Diag2(pixval,1)=Image(pixval,MdataSize+1-pixval);
% end
% Bckgrnd = [min(XAx),min(YAx),min(Diag1),min(Diag2)];
Z = Image;%-median(Bckgrnd);
% parameters are: [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]
x0 = Guess; % [Amp,x1,sigmax,y1,sigmay,0,Bckgrnd]; %Inital guess parameters
xs = x0; %centroid parameters
%noise = 0; % noise in % of centroid peak value (x(1))
InterpolationMethod = 'nearest'; % 'nearest','linear','spline','cubic'
FitForOrientation = 0; % 0: fit for orientation. 1: do not fit for orientation

%% ---Generate meshes for fitting--------------------------------------
xin = xs; 
% noise = noise/100 * x(1);
[X,Y] = meshgrid(x,y);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;
[Xhr,Yhr] = meshgrid(linspace(-MdataSize/2,MdataSize/2,300)); % generate high res grid for plot
xdatahr = zeros(300,300,2);
xdatahr(:,:,1) = Xhr;
xdatahr(:,:,2) = Yhr;

%% ---Gaussian Fit---------------------
if FitForOrientation == 0
    % define lower and upper bounds [Amp,xo,wx,yo,wy,fi]
%     lb = [                                  0,-MdataSize/4,           1, -MdataSize/4,           1, -pi/4];
%     ub = [realmax('double')-realmin('double'), MdataSize/4, MdataSize/2,  MdataSize/4, MdataSize/2,  pi/4];
    [xfit,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunctionRotConst,x0,xdata,Z,lb,ub);
else
    x0 =x0(1:5);
    xin(6) = 0; 
    xin =xin(1:5);
%     lb = [                                  0,-MdataSize/2,           1, -MdataSize/4,           1, -pi/4];
%     ub = [realmax('double')-realmin('double'), MdataSize/2, MdataSize/2,  MdataSize/2, MdataSize/2,  pi/4];
    [xfit,resnorm,residual,exitflag] = lsqcurvefit(@D2GaussFunction,x0,xdata,Z,lb,ub);
    xfit(6) = 0;
end

sizeParams=size(xfit);
GaussParams = zeros(sizeParams(1), sizeParams(2)+2); 
GaussParams(:,1:sizeParams(2)) = xfit;
Intensity = 0.9545*0.9545*2*3.14159*GaussParams(1).*GaussParams(3).*GaussParams(5); % This is the formula for the area within 2 sigma of the Gaussian;
GaussParams(:,sizeParams(2))=Intensity;
%Res = sum(sum(abs(residual)));%
Res = abs(sum(sum(residual)));% 
QuickPoisson = sqrt(mean(mean(Z))*3.14159*4.*GaussParams(3).*GaussParams(5)); % Estimate Poisson Noise from mean counts per image.
%QuickPoisson = Noise*3.14159*4.*GaussParams(3).*GaussParams(5);
GaussParams(:,sizeParams(2)+1)=Res+QuickPoisson; % Reasonable to question whether error bars should be 1 sigma or 2 sigma given noise. 

% Apply constraint. Amplitude must be greater than RoseCriterion*noise.
if GaussParams(:,1) <= RoseCriterion*Noise
    GaussParams(1) = 0;
    GaussParams(2) = NaN;
    GaussParams(4) = NaN;
    GaussParams(sizeParams(2))=0;
end

% R2 = RSquare2D(Image, residual)
% if R2 < 0.8
%    GaussParams(1)=0;
%    GaussParams(2)=NaN;
%    GaussParams(4)=NaN;    
% end

% GaussParams(4)

% %% ---------Plot 3D Image-------------
% figure(1)
% C = del2(Z);
% mesh(X,Y,Z,C) %plot data
% hold on
% surface(Xhr,Yhr,D2GaussFunctionRot(xfit,xdatahr),'EdgeColor','none') %plot fit
% axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5 ])%-noise noise+x(1)])
% alpha(0.2)  
% hold off
% 
% %% -----Plot profiles----------------
% hf2 = figure(2);
% set(hf2, 'Position', [20 20 950 900])
% alpha(0)
% subplot(4,4, [5,6,7,9,10,11,13,14,15])
% imagesc(X(1,:),Y(:,1)',Z)
% set(gca,'YDir','reverse')
% colormap('jet')
% 
% string1 = ['       Amplitude','    X-Coordinate', '    X-Width','    Y-Coordinate','    Y-Width','     Angle'];
% string2 = ['Set     ',num2str(xin(1), '% 100.3f'),'             ',num2str(xin(2), '% 100.3f'),'         ',num2str(xin(3), '% 100.3f'),'         ',num2str(xin(4), '% 100.3f'),'        ',num2str(xin(5), '% 100.3f'),'     ',num2str(xin(6), '% 100.3f')];
% string3 = ['Fit      ',num2str(xfit(1), '% 100.3f'),'             ',num2str(xfit(2), '% 100.3f'),'         ',num2str(xfit(3), '% 100.3f'),'         ',num2str(xfit(4), '% 100.3f'),'        ',num2str(xfit(5), '% 100.3f'),'     ',num2str(xfit(6), '% 100.3f')];
% 
% text(-MdataSize/2*0.9,+MdataSize/2*1.15,string1,'Color','red')
% text(-MdataSize/2*0.9,+MdataSize/2*1.2,string2,'Color','red')
% text(-MdataSize/2*0.9,+MdataSize/2*1.25,string3,'Color','red')
% 
% %% -----Calculate cross sections-------------
% % generate points along horizontal axis
% m = -tan(xfit(6));% Point slope formula
% b = (-m*xfit(2) + xfit(4));
% xvh = -MdataSize/2:MdataSize/2;
% yvh = xvh*m + b;
% hPoints = interp2(X,Y,Z,xvh,yvh,InterpolationMethod);
% % generate points along vertical axis
% mrot = -m;
% brot = (mrot*xfit(4) - xfit(2));
% yvv = -MdataSize/2:MdataSize/2;
% xvv = yvv*mrot - brot;
% vPoints = interp2(X,Y,Z,xvv,yvv,InterpolationMethod);
% 
% hold on % Indicate major and minor axis on plot
% 
% % % plot pints 
% % plot(xvh,yvh,'r.') 
% % plot(xvv,yvv,'g.')
% 
% % plot lins 
% plot([xvh(1) xvh(size(xvh))],[yvh(1) yvh(size(yvh))],'r') 
% plot([xvv(1) xvv(size(xvv))],[yvv(1) yvv(size(yvv))],'g') 
% 
% hold off
% axis([-MdataSize/2-0.5 MdataSize/2+0.5 -MdataSize/2-0.5 MdataSize/2+0.5])
% %%
% 
% ymin = 0;
% ymax = xfit(1);
% xdatafit = linspace(-MdataSize/2-0.5,MdataSize/2+0.5,300);
% hdatafit = xfit(1)*exp(-(xdatafit-xfit(2)).^2/(2*xfit(3)^2));
% vdatafit = xfit(1)*exp(-(xdatafit-xfit(4)).^2/(2*xfit(5)^2));
% subplot(4,4, [1:3])
% xposh = (xvh-xfit(2))/cos(xfit(6))+xfit(2);% correct for the longer diagonal if fi~=0
% plot(xposh,hPoints,'r.',xdatafit,hdatafit,'black')
% axis([-MdataSize/2-0.5 MdataSize/2+0.5 ymin*1.1 ymax*1.1])
% subplot(4,4,[8,12,16])
% xposv = (yvv-xfit(4))/cos(xfit(6))+xfit(4);% correct for the longer diagonal if fi~=0
% plot(vPoints,xposv,'g.',vdatafit,xdatafit,'black')
% axis([ymin*1.1 ymax*1.1 -MdataSize/2-0.5 MdataSize/2+0.5])
% set(gca,'YDir','reverse')
% %figure(gcf) % bring current figure to front

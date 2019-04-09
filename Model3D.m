

% This code assumes that you've already run the live script
% AtomAnalysisScript_Ethan_1298_trial and the peak positions and integrated
% intensities are located in the "RefinedPeaksGauss" variable.


    
% imageId = 105;  % image ID 
radius = 10;     % radius of 3D sphere
intensity = RefinedPeaksGauss(:,7);      % creates array of intensity from all values. Column 7 has integrated intensity data RefinedPeaksGauss(:,7).
intensityLimits = [min(intensity),max(intensity)]; % sets the limits as min/max of intensity
xLimits = [min(RefinedPeaksGauss(:,2))-10,max(RefinedPeaksGauss(:,2))+10];  
yLimits = [min(RefinedPeaksGauss(:,4))-10,max(RefinedPeaksGauss(:,4))+10];
fig = figure;

for imageId = 2 %select which image you want to look at for 3D model.
clf(fig);
[p] = find(RefinedPeaksGauss(:,10) == imageId);   %finds where imageID matches image # in the data

colormap jet
[xS, yS, zS] = sphere;


for pt = 1:numel(p)             % loops over each atom and plots it for each atom in each image 
%   for pt =15:15  
    n = p(pt);

    
    k = 15;
      for l = 100:100:RefinedPeaksGauss(n,7)
        
        hs=surf(RefinedPeaksGauss(n,2)+xS.*radius, ...  % Shift and scale x data
           RefinedPeaksGauss(n,4)+yS.*radius, ...  % Shift and scale y data
           k+zS.*radius, ...  % Shift and scale z data
           RefinedPeaksGauss(n,7).*ones(size(xS)), ... % color the sphere according to integrated intensity.
           'FaceColor', 'flat', ...
           'EdgeColor', 'none');
            k = k + 30;
            set(hs,'FaceColor','b')  % Uncomment this line to color each
                                         %  sphere a certain color.
       hold on;
      end  
   
end

% Modify figure and axes:
axis equal
xlim(xLimits);
ylim(yLimits);
zlim([0,120]);
set(gcf, 'Color', 'w');
%  set(gca, 'Color', 'w', 'Box', 'on', 'BoxStyle', 'back', ...
%       'XColor', 'none', 'YColor', 'none', 'ZColor', 'none', 'GridColor', 'none');
set(gca, 'ZTickMode', 'manual', 'Ztick', []);       % remove z axis ticks and labels
set(gca, 'XTickMode', 'manual', 'Xtick', []);       % remove x axis ticks and labels
set(gca, 'YTickMode', 'manual', 'Ytick', []);       % remove y axis ticks and labels

set(gca, 'Ydir', 'reverse');
colorbar;
colorbar('color','k');    % Change colorbar font to white
caxis(intensityLimits);         % set the color axis limits to be [min,max] of intensity values
% caxis([3000,12000]);
light('position', [2 -2 0.1], 'Style', 'infinite'); % create a light source
lighting phong;   % define how the light is reflected
%set(hs,'facealpha',0.8); % set the face alpha (opacity) to 80%
view([-3,71]);
% view([12,66]);
% save image as a .tif file. You need to copy and paste the path for
% saving.
% s = num2str(imageId);
% filename = strcat('CeO2cube6.3Dimage.',s, '.tif');
% fig.InvertHardcopy = 'off';
% figuresdir = 'C:\Users\ethan\Desktop\';
% saveas(fig, strcat(figuresdir,filename), 'tiff');

end
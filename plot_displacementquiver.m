function plot_displacementquiver(ZImage, P1, Displacements,scale)
%   plot_tracks lets you view an image series with the rmsds of atomic columns 
%   The rmsds are plotted on top of the Z-projection of an image series.
%   The rmsds are represented by markers on a colour scale, with different 
%   colours signifyign different magnitudes of displacement.  
%
%   INPUTS
%   
%   series - The image series in which you are tracking atom positions
%   
%   data - The output file from track_descrambler.m, containing lists of 
%   atomic column co-ordinates in different frames. 
%
%   Written by Barnaby Levin, ASU, 2017
sizeIM=size(ZImage);

x = P1(:,2);
y = P1(:,4);
u = Displacements(:,1);
v = Displacements(:,2);
ax = axes;
q = quiver(x,y,u,v,scale)
q.Color = 'red';
q.LineWidth = 1.5;
q.MaxHeadSize = 0.1;

set(ax, 'YDir', 'reverse'); % Ensure data plotted the right way up
daspect([1 1 1]); % Ensure data plotted with correct aspect ratio
xlim([0 sizeIM(2)]); % Ensure x axes of scatter plot and image are the same
ylim([0 sizeIM(1)]); % Ensure y axes of scatter plot and image are the same
set(gcf, 'color', [1 1 1]);
set(gca, 'xtick', []);
set(gca, 'ytick', []);

    % Comment/Uncomment this to see the ID numbers next to the peaks (this
    % will make diagnosis of bad peaks easier).
    b = num2str(round(Displacements(:,3),2))+" pm"; c = cellstr(b);
    dx = 3; dy = 3;
    text(P1(:,2)+dx, P1(:,4)+dy,c, 'Fontsize', 18, 'Color', 'k');
 

end
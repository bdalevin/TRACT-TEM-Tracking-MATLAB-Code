function plot_meanfwhm(ZProjImage, RefinedPeaks, Scale)
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

    sizeIM=size(ZProjImage);
    %Display Z-projected image
    ax1 = axes;
    imagesc(ax1, ZProjImage); axis image; axis off;
    frame_h = get(handle(gcf),'JavaFrame');
    set(frame_h,'Maximized',1);
    set(gcf,'color','w');
    
    ax2 = axes;
    % Plot the rmsds on top of the data for every atom 
    scatter(ax2, RefinedPeaks(:,2), RefinedPeaks(:,4), [], Scale*2355*(RefinedPeaks(:,3)+RefinedPeaks(:,5))/2, 'filled', 'LineWidth',1.5);
    % caxis ([30 60]);
    set(ax2, 'YDir', 'reverse'); % Ensure data plotted the right way up
    daspect([1 1 1]); % Ensure data plotted with correct aspect ratio
    xlim([0 sizeIM(2)]); % Ensure x axes of scatter plot and image are the same
    ylim([0 sizeIM(1)]); % Ensure y axes of scatter plot and image are the same
    linkaxes([ax1,ax2]) % Overlay the axes.
    ax2.Visible = 'off'; % Ensure image visible below scatter plot
    ax2.XTick = [];
    ax2.YTick = [];
    colormap(ax1,'gray') % Set colourmap for image to grayscale
    colormap(ax2,'jet') % Set colourmap for scatter plot to parula
    
    set([ax1,ax2],'Position',[.10 .11 .685 .815]); % Fix positions of the axes
    % These positions are in units of fractions of the canvas size. The
    % numbers are: Left of graph, Bottom of graph, width, and height. 
    % The width and heigh variables will be ineffective because of commands
    % from earlier in the code. 

    % Comment/Uncomment this to see the ID numbers next to the peaks (this
    % will make diagnosis of bad peaks easier).
    b = num2str(RefinedPeaks(:,6)); c = cellstr(b);
    dx = 3; dy = 3;
    text(RefinedPeaks(:,2)+dx, RefinedPeaks(:,4)+dy,c, 'Fontsize', 7, 'Color', 'w');
    
    cb2 = colorbar(ax2,'Position',[.78 .32 .0375 .4]); % Fix position of colourbar. 
    Ctitle = ylabel(cb2, 'Mean FWHM (pm)'); % Give colorbar a title

    % These positions are in units of fractions of the canvas size. The
    % numbers are: Left of colourbar, Bottom of colourbar, width, and height. 
    
end
function plot_gausspositions(ZImage, RefinedPeaks)
% Plot the lattice on top of the image so that the user can inspect it and
% make changes. 

% Make a figure to plot on
    h = figure('Name', 'Lattice Overlaid On Image', 'units','normalized','outerposition',[0 0 1 1]); 
    colormap('gray');

    %Display first image
    imagesc( ZImage ); axis image; 
    hold on; % Hold on forces every subsequent plot to be plotted on top of the image
    scatter(RefinedPeaks(:,2), RefinedPeaks(:,4),500,'x','r','linewidth', 2.0);  % Scatter plot the peaks found for the first image on top of the fist image
    brush on;
    
    % Comment/Uncomment this to see the ID numbers next to the peaks (this
    % will make diagnosis of bad peaks easier).
    b = num2str(RefinedPeaks(:,11)); c = cellstr(b);
    dx = 4; dy = 4;
    text(RefinedPeaks(:,2)+dx, RefinedPeaks(:,4)+dy,c, 'Fontsize', 16, 'Color', 'y');  
    
    
    hold off; % Allow the Figure to change again
end
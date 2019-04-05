function plot_gausspositions(ZImage, RefinedProjPeaksGauss)
% Plot the lattice on top of the image so that the user can inspect it and
% make changes. 

% Make a figure to plot on
    h = figure('Name', 'Lattice Overlaid On Image', 'units','normalized','outerposition',[0 0 1 1]); 
    colormap('gray');

    %Display first image
    imagesc( ZImage ); axis image; 
    hold on; % Hold on forces every subsequent plot to be plotted on top of the image
    scatter(RefinedProjPeaksGauss(:,2), RefinedProjPeaksGauss(:,4),'x','r','linewidth', 1.5);  % Scatter plot the peaks found for the first image on top of the fist image
    brush on;
    
%     % Comment/Uncomment this to see the ID numbers next to the peaks (this
%     % will make diagnosis of bad peaks easier).
%     b = num2str(RefinedProjPeaksGauss(:,6)); c = cellstr(b);
%     dx = 3; dy = 3;
%     text(RefinedProjPeaksGauss(:,2)+dx, RefinedProjPeaksGauss(:,4)+dy,c, 'Fontsize', 7, 'Color', 'w');
%   
    hold off; % Allow the Figure to change again
end
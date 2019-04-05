function plot_gausspositionsx2(ZImage, RefinedProjPeaksGauss1,RefinedProjPeaksGauss2)
% Plot the lattice on top of the image so that the user can inspect it and
% make changes. 

% Make a figure to plot on
    h = figure('Name', 'Lattice Overlaid On Image', 'units','normalized','outerposition',[0 0 1 1]); 
    colormap('gray');

    %Display first image
    imagesc( ZImage ); axis image; 
    hold on; % Hold on forces every subsequent plot to be plotted on top of the image
    scatter(RefinedProjPeaksGauss1(:,2), RefinedProjPeaksGauss1(:,4),'ro','filled','linewidth', 1.5);  % Scatter plot the peaks found for the first image on top of the fist image
    % Comment/Uncomment this to see the ID numbers next to the peaks (this
    % will make diagnosis of bad peaks easier).
    b = num2str(RefinedProjPeaksGauss1(:,6)); c = cellstr(b);
    dx = 3; dy = 3;
    text(RefinedProjPeaksGauss1(:,2)+dx, RefinedProjPeaksGauss1(:,4)+dy,c, 'Fontsize', 7, 'Color', 'r');
  
    
    scatter(RefinedProjPeaksGauss2(:,2), RefinedProjPeaksGauss2(:,4),'co','filled','linewidth', 1.5);
    % Comment/Uncomment this to see the ID numbers next to the peaks (this
    % will make diagnosis of bad peaks easier).
    b = num2str(RefinedProjPeaksGauss2(:,6)); c = cellstr(b);
    dx = -3; dy = 3;
    text(RefinedProjPeaksGauss2(:,2)+dx, RefinedProjPeaksGauss2(:,4)+dy,c, 'Fontsize', 7, 'Color', 'c');
    brush on;
  
    hold off; % Allow the Figure to change again
end
function plot_LatticevsGauss(ZImage, Lattice, RefinedPeaks)
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

    % Make a figure to plot on
    %h = figure('Name', 'Lattice Overlaid On Image', 'units','normalized','outerposition',[0 0 1 1]); 
    

    %Display first image
    imagesc( ZImage ); axis image; 
    colormap('gray');
    hold on; % Hold on forces every subsequent plot to be plotted on top of the image
    scatter(RefinedPeaks(:,2), RefinedPeaks(:,4), 'yo', 'filled', 'LineWidth',1.5);
    scatter(Lattice(:,1), Lattice(:,2),'x','b','linewidth', 1.5);  % Scatter plot the peaks found for the first image on top of the fist image

  
    hold off; % Allow the Figure to change again
    
end
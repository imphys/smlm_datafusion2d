%VISUALCLOUD2D Displays the binned density of the pointcloud in 2D.
% 
% SYNOPSIS:
%   visualizeCloud2D(pointcloud, bins, diameter, angle)
% 
% INPUT
%   pointcloud
%       [x y z; ...] (for example Particle{10}.points)
%   bins
%       number of bins in all directions
%   diameter
%       the diameter of the ROI in px
%   angle
%       rotate particle by this angle
% 
% OUTPUT
%   dip
%       dip image structure
%   Z
%       image matrix
%
% DEFAULTS:
%   none
% 
% NOTES:
% 
% (C) Copyright 2015               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Diederik Feilzer, July 2015
% Revision: Hamidreza Heydarian, 2017

function [dip, Z] = visualizeCloud2D(X, n, diameter, angle, caption)

    if nargin <5
        caption = '';
    end
    % rotation matrix
    Rot = [cos(angle) -sin(angle); sin(angle) cos(angle)];

    % apply rotation to particle
    X = (Rot*(X'))';

    % center particle around zero
    X(:,1) = X(:,1) - mean(X(:,1),1);
    X(:,2) = X(:,2) - mean(X(:,2),1);

    % extract the ROI
    ROIradius = 0.5*diameter;
    X = X(X(:,1) < ROIradius & X(:,1) > (-ROIradius),:);
    X = X(X(:,2) < ROIradius & X(:,2) > (-ROIradius),:);

    % define the grid
    xi = linspace(-ROIradius,ROIradius,n);
    yi = linspace(-ROIradius,ROIradius,n);

    % discretize the particle coordinates
    xr = interp1(xi,1:numel(xi),X(:,1),'nearest');
    yr = interp1(yi,1:numel(yi),X(:,2),'nearest');

    % binning
    Z = accumarray([xr yr],1, [n n]);
    
    % display the image with hot colormap
    dip = dipshow(Z','lin');
    dipmapping(dip,[0 max(Z(:))],'COLORMAP',hot(256));
    axis equal
    
    bar = 20*n/diameter/130; %20nm bar

    if caption
        text(n/2,n-10, caption,'Color','white','FontSize',15,'FontWeight' ...
             ,'bold','HorizontalAlignment','center','VerticalAlignment' ...
             ,'bottom');
    end
%     rectangle('Position',[n-140,n-70,bar,20],'FaceColor','white')    

end
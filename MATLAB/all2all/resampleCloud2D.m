% resampleCloud2D   subsamples the particle based on density of
% localizations
%
% SYNOPSIS:
%   [NewParticle] = resampleCloud2D(Particle)
%
% INPUT
%   Particle
%       The particle to be downsampled.
%
% NewParticle
%       The downsampled particle.
%
% NOTES
%       The weights for downsampling are computed using the intensity of 
%       the image resulting from binning and smoothing the localizations.
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Hamidreza Heydarian, 2017

function [ NewParticle ] = resampleCloud2D(Particle)

    meansigma = mean((Particle.sigma));         % the average uncertainties
    S = size(Particle.points,1);                % particles size
    cutoff = 5000;                              % the max number of 
                                                % localizations to be kept

    % binning the localizations
    xmax=max(Particle.points(:,1));
    xmin=min(Particle.points(:,1));
    ymax=max(Particle.points(:,2));
    ymin=min(Particle.points(:,2));
    dmax = [xmax ymax];
    dmin = [xmin ymin];
    nn = 100;                                   % number of bins
    bins = [1 1].*nn;
    binsize = (dmax - dmin)./bins;              % bin size

    xi = linspace(dmin(1),dmax(1),bins(1));
    yi = linspace(dmin(2),dmax(2),bins(2));

    xr = interp1(xi,1:numel(xi),Particle.points(:,1),'nearest');
    yr = interp1(yi,1:numel(yi),Particle.points(:,2),'nearest');

    subs = [xr yr];

    % act like addgaussianblob
    fsize = double(meansigma./binsize(1));  % filter size
    
    % smoothe the image
    f = fspecial('gaus',[round(7*fsize)+1 round(7*fsize)+1],fsize);
    
    % localizations is discretized to computed the weights for resampling
    binned = accumarray(subs,1, bins);          % discretized image
    localdensity = filter2(f,binned,'same');    % smoothed image

    % weights for resampling function
    weights = zeros(1,S);
    for l = 1:S
        weights(1,l) = localdensity(xr(l),yr(l));
    end
    
    % make sure that there is no NAN or negative weights
    weights(isnan(weights)) = 0;
    weights = max(weights,0);
    Max = numel(weights(weights > 0));

    % perform the weighted resampling
    ids = datasample(1:S,min(Max,cutoff),'Replace',false,'Weights',weights);

    % new particle
    NewParticle.points = Particle.points(ids,:);
    NewParticle.sigma = Particle.sigma(ids,:);

end
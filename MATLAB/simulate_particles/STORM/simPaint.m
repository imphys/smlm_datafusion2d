% simPaint   simulate SMLM particles from the model in simPaintModels
%
% SYNOPSIS:
%   [ coords, phantom, Q ] = simPaint( CCD_pixelsize, box_width, zoom, ...
%    timeframes, bg, n_photons, t_on, t_off, ...
%    model_name, raster_spacing, label_type, linker_size, dol, drift_factor,  ...
%    jit_binding_sites, verbose )
%
% INPUT
%     CCD_pixelsize     [nm per px]
%     box_width         [px]
%     zoom              [px per nm]
%     timeframes        Number of frames to measure [frames]
%     bg                background photons per pixel; TIRF 1-5, maybe 10, in 3D 50
%     n_photons         Expected collected photons per full frame
%     t_on              Mean lifespan of an on fluorescent label [frames]
%     t_off             Mean lifespan of an off fluorescent label [frames]
%     model_name        Model selection, see simPaintModels
%     raster_spacing    [nm]
%     label_type        Labeling type
%     linker_size       Size of fluorescent linkers [nm]
%     verbose           Verbose-mode
%
% OUTPUT
%   coords 
%       The structure that contains localization data
%
%   phantom
%       The phantom in image
%   binding_sites
%       The remainig binding sites
%
% NOTES
% 
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Author: Robert Nieuwenhuizen, Emile Boerder
% Edit:   Hamidreza Heydarian, 2017

function [coords, rigidParam, phantom, binding_sites, siteID] = simPaint( CCD_pixelsize, box_width, zoom, ...
   timeframes, bg, n_photons, t_on, t_off, ...
   model_name, raster_spacing, label_type, linker_size, dol, tr_mag, drift_factor,  ...
   jit_binding_sites, verbose, angle3D_factor )

%{
   Some great default parameters:
   
CCD_pixelsize = 160;           % [nm per px]
box_width = 512;               % [px]
zoom = 4;                      % [px per nm]

timeframes = 10000;            % Number of frames to measure [frames]
bg = 1;                        % background photons per pixel; TIRF 1-5, maybe 10, in 3D 50
n_photons = 1000;              % Expected collected photons per full frame
t_on = 3;                      % Mean lifespan of an on fluorescent label [frames]
t_off = 2e3;                   % Mean lifespan of an off fluorescent label [frames]

model_name = 'mpi';            % Model selection, see simPaintModels
raster_spacing = 5;            % [nm]
label_type = 'stick3D';        % Labeling type
linker_size = 3;               % Size of fluorescent linkers [nm]

verbose = 2;                   % Verbose-mode
   
%}

%% Just some parameters
if nargin<18
    angle3D_factor = 1;
end
fitboxsize = box_width/zoom/CCD_pixelsize;

lambda = 670;% emission color of marker
NA = 1.45;  
PSFsigma = 0.3*lambda/NA/CCD_pixelsize;

linker_size = linker_size / CCD_pixelsize;

init_Model = simPaintModels(model_name, raster_spacing, CCD_pixelsize);
if jit_binding_sites
    hypotenuse = 3/CCD_pixelsize;
    offsetX = (3/CCD_pixelsize)*rand(size(init_Model,1),1);
    offsetY = (hypotenuse^2 - offsetX.^2).^(0.5);
    init_Model = init_Model + [offsetX offsetY];
end
n_label_sites = size(init_Model,1);
Model = init_Model;

% first approach for removing binding sites
site_labeled = rand(n_label_sites,1)< dol;
Model(~site_labeled,:) = [];

% second approach for removing binding sites, assume gaussian distribution
% for the number of lablelled sites. The width of guassian is
% 0.1*number_of_binding_sites

% n_underlabel_sites = ceil(0.1*n_label_sites*randn + dol*n_label_sites); % \sigma*x+\mu
% site_idx = randperm(n_label_sites, n_label_sites - n_underlabel_sites);
% Model(site_idx,:) = [];

emitter_locs = Model - repmat(mean(init_Model),[size(Model,1) 1]);

%% investigating unbalanced rings
%{
idT = find(emitter_locs(:,3)>0);
idB = find(emitter_locs(:,3)<0);
rand_ang = (2*rand-1)*20;
emitter_locsT = emitter_locs(idT,:)*rotz(rand_ang)';
emitter_locs = [emitter_locsT;emitter_locs(idB,:)];
%}
%%
binding_sites = emitter_locs;   % remaining binding sites
%%

% determine the particle dimension 2D/3D
ndim = size(emitter_locs, 2);

if ndim == 2

    % filter outlier localizations (psuedo-segmentaiton)       
%         idx = find(coords(:,1) <  0.4  & ... 
%                   coords(:,1) > -0.4  & ...
%                   coords(:,2) <  0.4  & ...
%                   coords(:,2) > -0.4);    
%         coords = coords(idx,:);

%     idx = find(coords(:,1) <  1  & ... % for NUP107
%               coords(:,1) > -1  & ...
%               coords(:,2) <  1  & ...
%               coords(:,2) > -1);    
%     coords = coords(idx,:);

    particle_size = size(emitter_locs,1);


    % random translation
    translation = (tr_mag/CCD_pixelsize) * (rand(1,2)-0.5);
    translation_vec = repmat(translation, particle_size,1); 

    % random rotation
    angle = 2*pi*rand;
    rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)];

    % rotation + translation
    emitter_locs = emitter_locs * rot_matrix + translation_vec;
    
    % groundtruth parameters
    rigidParam = [angle translation];

else

    % filter outlier localizations (psuedo-segmentaiton)
%         idx = find(coords(:,1) <  0.65  & ...
%                   coords(:,1) > -0.65  & ...
%                   coords(:,2) <  0.65  & ...
%                   coords(:,2) > -0.65  & ...
%                   coords(:,3) <  0.40  & ...
%                   coords(:,3) > -0.40);    
%         coords = coords(idx,:);

%     idx = find(coords(:,1) <  0.8  & ...
%               coords(:,1) > -0.8  & ...
%               coords(:,2) <  0.8  & ...
%               coords(:,2) > -0.8  & ...
%               coords(:,3) <  0.6  & ...
%               coords(:,3) > -0.6);    
%     coords = coords(idx,:);

%     particle_size = size(emitter_locs, 1);

    % random translation
    translation = (tr_mag/CCD_pixelsize) * (rand(1,3)-0.5);

    % random rotation
%     angle = (2*rand(1,3)-1)*pi;             % [yaw, pitch, roll]
%     qtmp = angle2quat(angle(1), angle(2)/angle3D_factor, angle(3)/angle3D_factor);      % constraint yaw/pitch rotations % 2nd try +/-20 deg
    angle = genRandomAngle(angle3D_factor);    
    qtmp = angle2quat(angle(1), angle(2), angle(3));      % constraint yaw/pitch rotations % 2nd try +/-20 deg    
    q = [qtmp(2) qtmp(3) qtmp(4) qtmp(1)];

    % rotation + translation
    param = [q translation];

%     if i==1
%         param = [0 0 0 1 0 0 0];
%     end
    emitter_locs = transform_by_rigid3d(emitter_locs, param);

    % groundtruth parameters   
    rigidParam = param;

end

%% Switching simulation

% Simulate switching
N_emitters = size(emitter_locs,1);
M = poissrnd(timeframes/(t_on+t_off),N_emitters,1);

if exist('repelem', 'builtin')
    emitter_locs_on = repelem(emitter_locs, M, 1);
else
    emitter_locs_on = zeros(sum(M), size(emitter_locs, 2));
    j = 1;
    for i = 1:length(M)
        if M(i) == 0, continue, end
        emitter_locs_on(j:j+M(i)-1, :) = repmat(emitter_locs(i, :), M(i), 1);
        j = j + M(i);
    end
end

siteID = [];
for i=1:N_emitters
   siteID = [siteID;zeros(M(i),1)+i]; 
end

% Q = mean(M);

% Each frame has to be peturbed
frames_per_emitter = 1+poissrnd((t_on-1), length(emitter_locs_on), 1);

switch label_type
    
    case 'stick3D_projected'

        if size(emitter_locs_on, 2) > 2
            error(message('This model is for 2D data'));
        end
        
        % 3D stick with uniform distribution on sphere and the projection
        % on 2D surface
        for i = 1:length(emitter_locs_on)
            
            x = rand(frames_per_emitter(i),1);
            angles = acos(2*x-1);
            orientations = 2.*pi.*rand(frames_per_emitter(i),1);
            offset = linker_size * mean( ...
                [sin(angles) sin(angles)] .* ...
                    [sin(orientations) cos(orientations)]);
                
            emitter_locs_on(i, :) = emitter_locs_on(i, :) + offset;
        
        end

    case 'normal'
        
        dim = size(emitter_locs_on, 2);
        emitter_locs_on = emitter_locs_on + ...
            randn(size(emitter_locs_on)) .* ...
                repmat(linker_size ./ sqrt(frames_per_emitter), 1, dim);
            
    case 'stick3D'

        if size(emitter_locs_on, 2) < 3
            error(message('This model is for 3D data'));
        end
            
        x = rand(size(emitter_locs_on, 1), 1);
        phi = acos(2*x-1);
        theta = 2.*pi.*rand(1,1);

        xoff = cos(theta) .* sin(phi);
        yoff = sin(theta) .* sin(phi);
        zoff = cos(phi);
        denom = sqrt(xoff.^2 + yoff.^2 + zoff.^2);
        offset = linker_size * [xoff yoff zoff] ./ denom;
        emitter_locs_on = emitter_locs_on + offset;
        
end

%% Compute the CRLB on localization uncertainty

coords = emitter_locs_on;

% Crude approximation of number of photon distribution:
coords(:,4) = frames_per_emitter .* ...
    (1+geornd(1/n_photons, size(coords,1), 1));                                                     % Photon count
coords(:,8) = ceil(coords(:,4)/n_photons);                                                          % On frames
coords(:,6) = poissrnd(coords(:,8)*bg*fitboxsize^2)/fitboxsize^2;                                   % Background photons
coords(:,7) = PSFsigma*(1+0.02*randn(size(coords,1),1));                                            % Apparent PSF width
tau = 2*pi*(coords(:,7).^2+1/12).*coords(:,6)./coords(:,4);                                         % chemphyschem paper 2014 rieger & stallinga
coords(:,5) = sqrt((coords(:,7).^2+1/12)./coords(:,4).*(1+4*tau+sqrt((2*tau)./(1+4*tau))));         % Localization uncertainty, eq7 van ChemPhysChem
coords(:,9) = coords(:,5).*randn(size(coords(:,1)));
coords(:,1:2) = coords(:,1:2) + ...
                (coords(:,5)*[1 1]).*randn(size(coords(:,1:2))) + ...
                (drift_factor / CCD_pixelsize) * randn(size(coords(:,1:2)));                                                       % Add localization error

coords(:,10) = 3 * coords(:,5);                                                                     % localization uncertainty in z direction       
coords(:,3) = coords(:,3) + ...     % z direction
                (coords(:,10)).*randn(size(coords(:,3))) + ...
                (drift_factor / CCD_pixelsize) * randn(size(coords(:,3)));
            
% based on:  http://homepage.tudelft.nl/z63s8/Publications/Rieger2014.pdf

%remove errors
coords = coords(~any(isnan(coords)'),:);


%% For the verbose part

% Compute the phantom in image
phantom = [];   % default
if verbose > 1 && nargout > 2
    model_coords = floor( ( (Model + (fitboxsize/2)) ./ (fitboxsize) ) .* box_width );
    mean_coords = round(mean(model_coords) - (box_width / 2));
    model_coords = [model_coords(:, 1) - mean_coords(1), ...
        model_coords(:, 2) - mean_coords(2)];

    phantom = newim(box_width, box_width);
    for i = 1:length(model_coords)
        phantom(model_coords(i, 1), model_coords(i, 2)) = 1;
    end
end

% Send something to the screen
if verbose

    if verbose > 1

        if length(coords) < 5e4
            img = showlocalizations(coords(:,[1:2 5]) ./ fitboxsize, ...
                box_width, box_width, 'blobs', 3);
        else
            img = showlocalizations(coords(:,[1:2]) ./ fitboxsize, ...
                box_width, box_width);
        end
    
        img = addscalebar(img, 50, sprintf('%.1f nm', 50/zoom));
        fig = dipshow(img);
        colormap(fig, 'hot')
        hold on
        scatter(model_coords(:, 1), model_coords(:, 2), 'xc')
        hold off
        drawnow
        
    end
    
    fprintf('Binding sites: \t\t\t\t\t\t\t%d\n', length(Model))
    fprintf('Emitters per binding site: \t\t\t\t\t%.2f\n', N_emitters / length(Model))
    fprintf('Number of emitters: \t\t\t\t\t\t%.2e\n', N_emitters)
    fprintf('Number of localizations: \t\t\t\t\t%.2e\n', length(coords))
    fprintf('Q: \t\t\t\t\t\t\t\t%.2f\n', mean(M))
    fprintf('Raster spacing: \t\t\t\t\t\t%.2f nm\n', raster_spacing)
    fprintf('Label size: \t\t\t\t\t\t\t%.2f nm\n', linker_size*CCD_pixelsize)
    fprintf('Mean (median) localization uncertainty: \t\t\t%.2f nm\t(%.2f nm)\n', ...
        mean(coords(:,5))*CCD_pixelsize, median(coords(:,5))*CCD_pixelsize)
    fprintf('Mean number of photons (after frame conn): \t\t\t%.2e\n', mean(coords(:, 4)))
    fprintf('Mean number of bg photons (per pixel after frame conn): \t%.2e\n', mean(coords(:, 6)))

end
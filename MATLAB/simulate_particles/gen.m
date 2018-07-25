% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Author: Emile Boerder, 2016
% Edit:   Hamidreza Heydarian, 2017
    
%% The number of particles
N = 256;

%% Configuration of simPaint
CCD_pixelsize = 130;                    % [nm per px]
box_width = 512;                        % [px]
zoom = 4;                               % [px per nm]

timeframes = 100000;                    % Number of frames to measure [frames]
bg = 1;                                 % background photons per pixel; TIRF 1-5, maybe 10, in 3D 50
n_photons = 5000;                       % Expected collected photons per full frame
dol = 0.7;                              % density of labeling 
t_on = 3;                               % Mean lifespan of an on fluorescent label [frames]
t_off = 2e3;                            % Mean lifespan of an off fluorescent label [frames]

model_name = 'tud_ralf';                % Model selection, see simPaintModels
raster_spacing = 5;                     % [nm]
label_type = 'stick3D';                 % Labeling type
linker_size = 0.66;                     % Size of fluorescent linkers [nm]
drift_factor = 0;                       % drift factor [nm]
jit_binding_sites = 0;                  % flag for jittering the binding site to have then out of grid
tr_mag = 10;                            % translation magnitude [nm]       

verbose = 0;                            % Verbose-mode
filename = ['tud_nph' ...
            num2str(n_photons) ...
            '_dol' ...
            num2str(uint8(100*dol)) ...
            '_tr' ...
            num2str(tr_mag) ...
            '_nm' ...
            'N' ...
            num2str(N) ...
            'sym'];                     % file name to save

%% Creation of the particles
Particles = cell(1, N);

for i = 1:N

   coords = simPaint( CCD_pixelsize, box_width, zoom, ...
       timeframes, bg, n_photons, t_on, t_off, ...
       model_name, raster_spacing, label_type, linker_size, dol, drift_factor, ...
       jit_binding_sites, verbose );
   
   idx = find(coords(:,1) <  0.4  & ...
              coords(:,1) > -0.4  & ...
              coords(:,2) <  0.4  & ...
              coords(:,2) > -0.4);    
   coords = coords(idx,:);
   particle_size = size(coords(:,1),1);
   translation = (tr_mag/CCD_pixelsize) * (rand(1,2)-0.5);
   translation_vec = repmat(translation,particle_size,1); 
   angle = 2*pi*rand;
   rot_matrix = [cos(angle) -sin(angle); sin(angle) cos(angle)];
   coords(:, 1:2) = coords(:, 1:2) * rot_matrix + translation_vec;
    
   Particles{1,i} = struct('coords', coords, 'angle', angle, 'translation', translation);
    
end


%% Visualize one random particle for inspection
% lucky = ceil(N * rand);
% visualizeCloud2D(Particles{1,lucky}.coords(:,1:2), 600, 0.6, 0, 'Sample simulated particle');

% scatter plot
% figure;
% scatter(Particles{1,lucky}.coords(:,1),Particles{1,lucky}.coords(:,2),'.')
% axis square

% 
% fitboxsize = box_width/zoom/CCD_pixelsize;
% img = showlocalizations(...
%     Particles{lucky}.coords(:,[1:2 5]) ./ fitboxsize, box_width, box_width, ...
%         'blobs', 3);
%     
% img = addscalebar(img, 50, sprintf('%.1f nm', 50/zoom));
% fig = dipshow(img);
% colormap(fig, 'hot')


%% Save everything
if ~exist('data/test','dir')
    mkdir('data/test');
end
save(['data/test/' filename], 'Particles');

%% some analysis
% nlocs = zeros(1,N);
% for i=1:N
%     nlocs(i) = size(Particles{1,i}.coords(:,1:2),1);
% end
% meannlocs = mean(nlocs)

%% compute the ground truth for simulated data set (TUD)

% superParticle = [];
% N = numel(Particles);
% for i = 1:N
%     
%     particle_size = size(Particles{1,i}.coords(:,1),1);
%     ang = Particles{1,i}.angle;
%     tr = Particles{1,i}.translation;
%     points = Particles{1,i}.coords(:,1:2);
%     R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
%     
%     tmpParticle = (points - repmat(tr,particle_size,1))*R';
%     superParticle = [superParticle; tmpParticle];
%     
% end
% 
% visualizeCloud2D(superParticle,600,0.6,0,'Ground-truth');
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
clc
close all
N = 1;

%% Configuration of simPaint
CCD_pixelsize = 130;                    % [nm per px]
box_width = 512;                        % [px]
zoom = 4;                               % [px per nm]

dol = 1;                                % density of labeling 

% model_name = '2D_NUP107';               % Model selection, see simPaintModels
model_name = '2D_ringTeun'; 
raster_spacing = 5;                     % [nm]
label_type = 'stick3D_projected';                 % Labeling type. 2D: 'stick3D_projected', 3D: 'stick3D'
linker_size = 2*0.66;                     % Size of fluorescent linkers [nm]
drift_factor = 0;                       % drift factor [nm]
jit_binding_sites = 0;                  % flag for jittering the binding site to have them out of grid
tr_mag = 20;                            % translation magnitude [nm]   
ndim = 2;                               % 2D or 3D data (=2 or 3)
angle3D_factor = 3;                   % the factor to limit rotation around x,y

% palm, storm, gsdim settings
% Parameters related to switching statistics
switching_params = struct();

% switching_params.tech = 'paint';
% switching_params.timeframes = 0.5E5;      % paint
% t_on = 3;                               % Mean lifespan of an on fluorescent label [frames] storm: 3
% t_off = 2e3;                            % Mean lifespan of an off fluorescent label [frames] storm: 2e2

switching_params.tech = 'storm';
switching_params.timeframes = 1E3;      % storm:1E3
t_on = 3;                               % Mean lifespan of an on fluorescent label [frames] storm: 3
t_off = 2e2;                            % Mean lifespan of an off fluorescent label [frames] storm: 2e2

switching_params.k_on = 1/t_on;         %1E-3;
switching_params.k_off = 1/t_off;       %1;
switching_params.k_bleach = 0.02;
switching_params.subframes = 5;

% Parameters related to localization statistics
localization_params = struct();
localization_params.photons = 1000;      % Expected collected photons per full frame (photon count)
localization_params.bg = 300;             % background photons per pixel; TIRF 1-5, maybe 10, in 3D 50
localization_params.PSFsigma_xy = 580/4/1.49;

verbose = 0;                            % Verbose-mode
%%
row_names = {'N=';'CCD_pixelsize=';'box_width=';'zoom=';'dol=';'t_on=';'t_off='; ...
    'raster_spacing=';'linker_size=';'drift_factor='; ...
    'jit_binding_sites=';'tr_mag=';'ndim='; ...
    'switching_params.timeframes='; 'switching_params.k_on='; ...
    'switching_params.k_off='; 'switching_params.k_bleach='; ...
    'switching_params.subframes='; 'localization_params.photons='; ...
    'localization_params.bg='; 'localization_params.PSFsigma_xy='; ...
    'angle3D_factor'};
T = table([ ...
    N;CCD_pixelsize;box_width;zoom;dol;t_on;t_off; ...
    raster_spacing;linker_size;drift_factor; ...
    jit_binding_sites;tr_mag;ndim; ...
    switching_params.timeframes;switching_params.k_on; ...
    switching_params.k_off;switching_params.k_bleach; ...
    switching_params.subframes;localization_params.photons; ...
    localization_params.bg;localization_params.PSFsigma_xy;angle3D_factor],'RowNames',row_names);

if switching_params.tech == 'storm'
    
    filename = [num2str(N) '_' model_name ...
                '_nph' ...
                num2str(localization_params.photons) ...
                '_dol' ...
                num2str(uint8(100*dol)) ...
                '_' ...
                num2str(ndim) ...
                'D_' ...
                switching_params.tech ...
                '_ang' ...
                num2str(floor(10*angle3D_factor))];                     % file name to save    
else % storm
    
    filename = [num2str(N) '_' model_name ...
                '_nph' ...
                num2str(localization_params.photons) ...
                '_dol' ...
                num2str(uint8(100*dol)) ...
                '_' ...
                num2str(ndim) ...
                'D_' ...
                switching_params.tech ...
                '_ang' ...
                num2str(floor(10*angle3D_factor))];                     % file name to save    
    
end

%% Creation of the particles
particles = cell(1, N);

for i = 1:N

    if switching_params.tech == 'paint'
       [coords, rigidParam, ~, binding_sites, siteID] = simPaint( CCD_pixelsize, box_width, zoom, ...
           switching_params.timeframes, localization_params.bg, localization_params.photons, t_on, t_off, ...
           model_name, raster_spacing, label_type, linker_size, dol, tr_mag, drift_factor, ...
           jit_binding_sites, verbose, angle3D_factor);

        if ndim == 3
            
            idx = find(coords(:,1) <  0.8  & ...
                      coords(:,1) > -0.8  & ...
                      coords(:,2) <  0.8  & ...
                      coords(:,2) > -0.8  & ...
                      coords(:,3) <  0.6  & ...
                      coords(:,3) > -0.6);    
            coords = coords(idx,:); 
            siteID = siteID(idx,:);
            param = rigidParam;
            particles{1,i} = struct('coords', coords, 'param', param, 'binding_sites', binding_sites, 'binding_sites_ID', siteID); 
            
        else
            
            idx = find(coords(:,1) <  1  & ... % for NUP107 1; for tud 0.4
                       coords(:,1) > -1  & ...
                       coords(:,2) <  1  & ...
                       coords(:,2) > -1);    
            coords = coords(idx,:);
            angle = rigidParam(1);
            translation = rigidParam(2:3);

            particles{1,i} = struct('coords', coords, 'angle', angle, ...
                                   'translation', translation);   
                               
        end      
       
    else
       [coords, rigidParam, ~, binding_sites, siteID] = simPalmStorm( CCD_pixelsize, box_width, zoom, ...
           switching_params.timeframes, localization_params.bg, localization_params.photons, t_on, t_off, ...
           model_name, raster_spacing, label_type, linker_size, dol, tr_mag, drift_factor, ...
           jit_binding_sites, verbose, switching_params, localization_params, angle3D_factor); 
       
        if ndim == 3
            
            idx = find(coords(:,1) <  0.8  & ...
                      coords(:,1) > -0.8  & ...
                      coords(:,2) <  0.8  & ...
                      coords(:,2) > -0.8  & ...
                      coords(:,3) <  0.6  & ...
                      coords(:,3) > -0.6);    
            coords = coords(idx,:);  
            siteID = siteID(idx,:);
            param = rigidParam;
            particles{1,i} = struct('coords', coords, 'param', param, 'binding_sites', binding_sites, 'binding_sites_ID', siteID); 
            
        else
            
            idx = find(coords(:,1) <  1  & ... % for NUP107 1; for tud 0.4
                       coords(:,1) > -1  & ...
                       coords(:,2) <  1  & ...
                       coords(:,2) > -1);    
            coords = coords(idx,:);
            angle = rigidParam(1);
            translation = rigidParam(2:3);

            particles{1,i} = struct('coords', coords, 'angle', angle, ...
                                   'translation', translation);   
                               
        end            

    end 
    
end

% exclude null particles
for i=1:N
    sizen(i) = size(particles{1,i}.coords,1);
end
particles(sizen==0) = [];
%% Visualize one random particle for inspection

lucky = ceil(N * rand);

if ndim == 2 
    visualizeCloud2D(particles{1,lucky}.coords(:,1:2), 100, 2.5, 0, ...
                    'Sample simulated particle');
else
    figure;
    scatter3(particles{1,lucky}.coords(:,1),particles{1,lucky}.coords(:,2) ...
            ,particles{1,lucky}.coords(:,3),'.');
        axis square;axis equal;
end

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
% curDir = cd;
% if ~exist([curDir '/data'],'dir')
%     mkdir('data');
% end
% % 
% mkdir(['data/' filename])
% save(['data/' filename '/' filename], 'particles');
% writetable(T,['data/' filename '/' filename '.txt'],'Delimiter',' ','WriteRowNames',true);

%% some analysis
N = numel(particles);
nlocs = zeros(1,N);
for i=1:N
    nlocs(i) = size(particles{1,i}.coords(:,1:2),1);
end
meannlocs = mean(nlocs);
% locspersite = meannlocs / 37      %depends on DoL!
for i=1:N
    sigmas(i) = mean(particles{1,i}.coords(:,5));
    
end

disp(['average number of locs per particle is ' num2str(meannlocs)]);
disp(['mean localization unc is ' num2str(CCD_pixelsize*mean(sigmas))]);

%% compute the ground truth for simulated data set (TUD)

% ndim = 3;
% CCD_pixelsize = 130; 

N = numel(particles);
if ndim == 2
%     superParticle = [];
%     
%     for i = 1:N
% 
%         particle_size = size(particles{1,i}.coords(:,1),1);
%         ang = particles{1,i}.angle;
%         tr = particles{1,i}.translation;
%         points = particles{1,i}.coords(:,1:2);
%         R = [cos(ang) -sin(ang); sin(ang) cos(ang)];
% 
%         tmpParticle = (points - repmat(tr,particle_size,1))*R';
%         superParticle = [superParticle; tmpParticle];
% 
%     end
% 
%     visualizeCloud2D(superParticle,300,2,0,'Ground-truth');
    
else % visualize super 3D particle
    
    sr = [];
    for i=1:numel(particles)

        % fetch the transformation parameters
        tmpParam = particles{1,i}.param;

%         % first translate
        gtParam1 = [0, 0, 0, 1, ...
                   -tmpParam(5),-tmpParam(6),-tmpParam(7)];           
        tmpParticle = transform_by_rigid3d(particles{1,i}.coords(:,1:3), gtParam1);

        % second rotate
        gtParam2 = [-tmpParam(1),-tmpParam(2),-tmpParam(3),tmpParam(4), 0, 0, 0];
        tmpParticle = transform_by_rigid3d(tmpParticle, gtParam2);

        sr = [sr;tmpParticle];

    end

    % visualize
    sr = CCD_pixelsize*sr;
%     figure
%     scatter3(sr(:,1),sr(:,2),sr(:,3),'.')
    visualizeSMLM3D(sr,5);
    title(['Ground-truth'])
    
    sup = [];
    for i=1:N
        sup=[sup;particles{1,i}.coords(:,1:3)];
    end    
    visualizeSMLM3D(sup,0.025);
    title(['Accumulation'])    
end


%%
% s = 0;
% for i=1:size(M,1)
%     if i==1
%         t(1:M(i),2) = i;
%     else
%         s = s + M(i-1);
%         t(s+1:s+M(i),2) = i;
%     end
% end
% %%
% tmax = max(t(:,1));
% for i=1:10
%     idx = find(t(:,1)>(i-1)*100 & t(:,1)<i*100);
%     tmp{1,i} = t(idx,2);
% end
% %%
% for i=1:10
%     for j=1:37
%         id = find(tmp{1,i}==j);
%         hh(i,j) = size(id,1);
%     end
% end
% Define Markov model
% k_on
% 
% trans_mat = [(1-k_on/subframes) k_on/subframes 0; k_off/subframes 1-(k_off+k_b)/subframes k_b/subframes; 0 0 1];
% visible_mat = [1 0;0 1;1 0];
% 
% % Simulate switching
% emitter_states = [];
% 
% for nn = 1:N_emitters
%     [seq states] = hmmgenerate(timeframes*subframes,trans_mat,visible_mat);
%     seq = seq - 1;
%     seq = sum(reshape(seq,[subframes timeframes]),1);
%     % Format of emitter_states: [emitter_id x y t n_photons]
%     if (sum(seq~=0) ~= 0)
%         emitter_states = cat(1,emitter_states,[ones(sum(seq~=0),1)*[nn emitter_locs(nn,:)] find(seq)' nonzeros(seq)]);
%     end
% end
% This script loads a dataset from '/data' directory which should be a v
% cell array of structures called particles with localization in 'point' 
% field (x,y) and uncertainties in 'sigma' field:
%   particles{1,1}.points        -> localization data (x,y) in camera pixel
%   particles{1,1}.sigma         -> localization uncertainties (sigma) in camera pixel
%
% The experimental data which are stored in 'data/experiment/' directory have 
% already this format.
%
% The simulated data which are stored in 'data/simulation/' have a
% different format. Each single simulated particle has the following format:
%   particles{1,1}.coords(:,1:2) -> localization data (x,y) in camera pixel
%   particles{1,1}.coords(:,5)   -> localization uncertainties (sigma)
%   particles{1,1}.angle         -> ground truth angle in radian
%   particles{1,1}.translation   -> ground truth translation vector (tx,ty)
%   in camera pixel.
% Please refer to gen.m in 'MATLAB/simulate_particles/' for the descriptions 
% of other fields in particles{1,1}.coords.
%
% The following code will load (experimental/simulated) data from file and
% then performs the 3 steps particle fusion algorithm. Different datasets
% are provided. You only need to provide the 'dataset' name and the number 
% of particles to be averaged (K).
%
% The code for simulating SMLM particles is also provided in the
% 'MATLAB\simulate_particle' directory. You need to run the gen.m script to
% make simulated particles. The generated particles are
% stored in the 'data/test' directory which can again be used with the current
% script.
%
% The codes makes use of parralel computing toolbox to distribute the load
% over different workers. You increase this above the fixed value (12) by MATLB
% by changing your parallel.settings.txt file in the ./matlab/verision/
% directory.
% <key name="PreferredNumWorkers"><double><value>48.0</value></double>
% The code will automatically use the GPU if available.


% (C) Copyright 2018                    Deparment of Imaging Physics
%     All rights reserved               Faculty of Applied Sciences
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
% Hamidreza Heydarian <H.Heydarian@tudelft.nl>
% Bernd Rieger <B.Rieger@tudelft.nl>

%%
close all
clear all
%add the required directory to path
path_matlab = genpath('MATLAB');
addpath(path_matlab)

%% LOAD DATASET

% load dataset stored in data directory
% datasets used in the paper 
% simulation | experiment 
% tud_dol100 | tud_dol80_raw
% tud_dol50  | tud_dol80_filtered
% tud_dol30  | tud_dol50_raw
%            | tud_dol50_filtered
%            | tud_dol30

% Set the dataset parameters
% --- experimental data sets ---
%dataset = 'experiment/tud_dol80_filtered';          
%dataset = 'experiment/tud_dol50_filtered';          
%dataset = 'experiment/tud_dol80_raw';          
%dataset = 'experiment/tud_dol50_raw';

% ---- simulated data sets ---
 dataset = 'simulation/tud_dol100.mat';             
% dataset = 'simulation/tud_dol50.mat';             
% dataset = 'simulation/tud_dol30.mat';  

K = 50;     % set to 1e3 for all particless        % choose K particles

% load the dataset
load(['data/' dataset])

% decide whether the dataset is from experiment or from simulation
flag = isfield(particles{1,1}, 'points');   % 1/0 for experiment/simulation

if K > numel(particles)
    K = numel(particles);
end

subParticles = cell(1,K);

if flag
    % experimental dataset
    for j=1:K
        subParticles{j}.points = double(particles{1,j}.points);
        subParticles{j}.sigma  = double((particles{1,j}.sigma).^2);

    end
else    
    % simulated dataset
    for j=1:K
       subParticles{j}.points = particles{1,j}.coords(:,1:2);
       subParticles{j}.sigma  = particles{1,j}.coords(:,5).^2;
    end
end

%% STEP 1

% perform all2all registration and save results on disk
scale = 0.01;   % For TUD logo, a value of 0.01 (in camera pixel unit) 
                % works well. For NPC data, 0.1 is giving better result.
                % Look at Online Methods for the description
all2all(subParticles, 'output/all2all_matrix', scale);

%% STEP 1,2

% Lie-algebraic averaging and outlier removals (step 1,2)
[initAlignedParticles, M1] = outlier_removal(subParticles, 'output/all2all_matrix/', 'output/');

%% STEP 3

% bootstrapping (step 3)
iter = 4;   % number of iterations
% [superParticle, M2] = one2all(initAlignedParticles, iter, M1, 'output/');  

% for NPC data use the following function which takes into account the
% prior knowledge about 8-fold symmetry. Look at Online Methods for details
[superParticle, M2] = one2all_symmetric(initAlignedParticles, iter, M1, 'output/');

%% VISUALIZE RRESULTS

% visualize a random sample particle for inspection
lucky = ceil(K * rand);
visualizeCloud2D(subParticles{1,lucky}.points, 600, 0.6, 0, 'Sample initial particle');

% visualize the particle fusion reconstruction at iter=2. Change to >=2 for 
% reconstruction at further iterations
iter = 2;                                            
visualizeCloud2D(superParticle{1,iter},600,0.6,0, 'Reconstruction');

% visualize grount-truth for the subset of K particles for simulated data
if ~flag
    groundtruth = [];
    for i = 1:K

        curSize = size(particles{1,i}.coords(:,1),1);   % current particle size
        ang = particles{1,i}.angle;                     % ground-truth angles
        tr = particles{1,i}.translation;                % ground-truth translation
        points = particles{1,i}.coords(:,1:2);          % localizations
        R = [cos(ang) -sin(ang); sin(ang) cos(ang)];    % rotation matrix

        tmpParticle = (points - repmat(tr,curSize,1)) * R';
        groundtruth = [groundtruth; tmpParticle];       % ground-truth reconstruction

    end
    visualizeCloud2D(groundtruth,600,0.6,0,'Ground-truth');
end
% pairFitting_symmetric   register a pair of particles (parallelized with
% parfor). At the end, the particle is rotated with a random multiple of 
% 2*pi/8 angle to distribute evenly on the 8-blobs. Refer to the paper (figure 3)
%
% SYNOPSIS:
%   [parameter, registered_model, history, config, max_value] = pairFitting_symmetric(M, S, scale)
%
% INPUT
%   M
%       The first particle.
%   S
%       The second particle.nano 
%   scale
%       The array containing the set of scale. The size of the array should
%       be equal to iter in one2all_symmetric() function
%
% OUTPUT
%   parameter
%       Rigid registration parameter [angle, t1, t2].
%   registered_model
%       The result of registering M to S.
%   history
%       History of the optimization variables over iterations (not used!).
%   config
%       The structure containing all the input parameters for GMM based 
%       registration, see '/MATLAB/initialize_config.m'.
%   max_value
%       The maximum value of the Bhattacharya cost function.
%
% NOTES
%       Registration algorithms normally work for a certain range of 
%       rotation angle and scales. In order to avoid trapping in local 
%       minima, different initializations is provided. Then, GMM-based 
%       method, registers the two particles using different
%       intializations. Finally, Bhattacharya cost function, chooses the
%       rigid parameter set which gives the highest score.    
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Hamidreza Heydarian, 2017

function [parameter, registered_model, history, config, max_value] = pairFitting_symmetric(M, S, scale) 
                                    
    Nfold = 8; 
    angle = linspace(-pi,pi-2*pi/Nfold,Nfold);  %initial angles
    
	% generate all permutation of initial grid and scale
	[init_scale, init_angle] = ndgrid(scale, angle);
	N_init = numel(init_scale);

    % initalize functions outputs
    param = cell(1, N_init);                    % rotation & translation
    transformed_model = cell(1, N_init);        % fused particle
    bhattacharyya_distance = zeros(1, N_init);  % bhattacharyya distance
    history = cell(1, N_init);                  % optimization variable
    config = cell(1, N_init);                   % struct for GMM algorithm
    
    % resample particles to save time
    M_resampled = resampleCloud2D(M);
    %S_resampled = resampleCloud2D(S); 
    
    %during Bootstrapping: resample S randomly, otherwise bright spots are
    %enhanced
    iS = randsample(size(S.points,1),min([5000,size(S.points,1)])); 
    S_resampled.points = S.points(iS,:);
    S_resampled.sigma = S.sigma(iS); 

	% loop over scales and inital angles
    parfor iter = 1:N_init
        
        % initialize the GMM registration method
        f_config = initialize_config(M_resampled.points, S_resampled.points, 'rigid2d');
        f_config.init_param = [0, 0, init_angle(iter)];
        f_config.scale = init_scale(iter);

        % perform registration
        [param{iter}, transformed_model{1, iter}.points, history, config, fval(iter)] = gmmreg_L2(f_config);
        transformed_model{1, iter}.sigma = [M_resampled.sigma; S_resampled.sigma];
        
        bhattacharyya_distance(iter) = expdist(S_resampled, M_resampled, param{iter}(1,3), param{iter}(1,1), param{iter}(1,2));
        
    end
			
    % maximize over scales and inital angles
    [max_value, IDX] = max(bhattacharyya_distance);

    % random rotation by integer multiple of 2*pi/8
    finalAngleID = randperm(size(fval,2),1);
    parameter = param{1,IDX};
    parameter = [0, 0, angle(finalAngleID)+parameter(1,3)];

    M_transformed.points = transform_by_rigid2d(M.points, parameter);
    M_transformed.sigma = M.sigma;
    registered_model.sigma = [M_transformed.sigma; S.sigma]; 
    registered_model.points = [M_transformed.points; S.points]; 
	        
end

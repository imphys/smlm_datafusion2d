%outlier_removal   This function performs Lie-algebraic averaging and
% registration outlier removal and then performs the second round of
% Lie-algebraic averaging
%
%   SYNOPSIS:
%       [Particles_step2, M_new] = outlier_removal(particles, all2all_dir, outdir)
%
%   Input: 
%       particles: Cell array of particles of size 1xN
%       all2all_dir: The directory in which rows of all2all registration
%       matrix is stored
%       outdir: Output directory where the results are stored
%
%   Output:
%       Particles_step2: the aligned particles after outlier removal
%       M_new: Transformation parameters (rotation+translation), a 4x4xN
%       matrix where N is the number of particles.
%
%   NOTE:
%
%
% (C) Copyright 2017                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Hamidreza Heydarian, Oct 2017.

function [initAlignedParticles, M_new] = outlier_removal(particles, all2all_dir, outdir)

    disp('Lie-algebraic averaging started  !');
    path_matlab = genpath('Optimization');
    addpath(path_matlab)
%     cvx_solver mosek;
%     cvx_solver

    % initialization
    RM = zeros(4,4,5);
    I = zeros(2,5);
    iter = 1;   

    % load all2all registration matrix rows from file
    allRows = dir(all2all_dir);
    allNames = {allRows(~[allRows.isdir]).name};
    allNames = natsortfiles(allNames);
    nRows = numel(allNames);

    % stack all relative motion parameteres in RM
    for j=1:nRows

        load([all2all_dir allNames{1,j}]);
        N = numel(result);

        for k=1:N

            I(:,iter) = result(k).id;

            estAngleR = -result(k).parameter(3);
            curestAngle = wrapToPi(estAngleR);
            w_est = curestAngle * [0;0;1];
            RM(1:3,1:3,iter) = expm([0,-w_est(3),w_est(2); w_est(3),0,-w_est(1); -w_est(2),w_est(1),0]);
            RM(:,4,iter) = [result(k).parameter(1:2)'; 0; 1];                

            iter = iter + 1;
        end

    end

    nParticles = numel(particles);

    % perform Lie-algebraic averagin
    totalPairs = 0.5 * nParticles * (nParticles-1);
    edges = motion_graph(nParticles);
    for i=1:nParticles-1
        idxPair(i) = find(I(1,:) == edges(1,i) & I(2,:) == edges(2,i));
    end
    nPairs = totalPairs - numel(idxPair); % change to average only a subset
    totalPairsIDX = 1:totalPairs;
    totalPairsIDX(idxPair) = [];
    idxRest = randperm(numel(totalPairsIDX), nPairs);
    idxPair = [idxPair totalPairsIDX(idxRest)];
    [M] = MeanSE3Graph(RM(:,:,idxPair),I(:,idxPair));

    % the aligned particles from Lie-algebraic averaging (step 1)
    Particles_step1 = applyRigidTransform(particles, M); 

    % initial relative all2all rotations    
    relTr_init = zeros(2,5);
    for i=1:totalPairs
        relAngleR_init(i) = atan2(RM(2,1,i),RM(1,1,i));
        relTr_init(:,i) = RM(1:2,4,i);
    end

    % relative all2all rotation and translation after averaging
    kk = 1;
    relTr = zeros(2,5);
    for i=1:nParticles-1
        for j=i+1:nParticles

            relAngleR(kk) = atan2(M(2,1,i),M(1,1,i)) - atan2(M(2,1,j),M(1,1,j));
            relTr(:,kk) = M(1:2,4,i)-M(1:2,4,j);
            kk = kk+1;

        end
    end
    
    disp('Outlier removal and second Lie-algebraic averaging started  !');
    % registration outlier removal
    error = wrapTo180(rad2deg(relAngleR + relAngleR_init));
    error_tr = sqrt(sum((relTr - relTr_init).^2,1));    % translation outliers, 
                                                        % not used
    error_idx = find(abs(error) > 5);                   % rotation outliers, the
                                                        % threshold is set
                                                        % to 5 degree
    
    RM_new = RM;
    I_new = I;
    RM_new(:,:,error_idx) = [];                         % exclude outliers                        
    I_new(:,error_idx) = [];                            % exclude outliers

    if isConnected(I_new)==0    
        
        error_idx = connectGraph(I,error);
        
        RM_new = RM;
        I_new = I;
        RM_new(:,:,error_idx) = [];                         % exclude outliers                        
        I_new(:,error_idx) = [];  
        disp('Reconnected Graph')
    end
    
    % second round of Lie-averaging
    [M_new] = MeanSE3Graph(RM_new,I_new);

    % the superparticle from Lie-algebraic averaging after outlier removal (step 2)
    initAlignedParticles = applyRigidTransform(particles, M_new); 

    % save lie-averaging and outlier removal output to file
    save([outdir '/motion_Lie_averaging'], 'M');
    save([outdir '/particle_Lie_averaging'], 'Particles_step1');

    save([outdir '/motion_outlier_removal_'], 'M_new');
    save([outdir '/particle_outlier_removal_'], 'initAlignedParticles');
    
    disp('Lie-algebraic averaging and outlier removal is done !');
    fprintf('\n\n');
    
end
%one2all   Register each particle to the stack of all the remaining
%particles and iterate
% of particles
%
%   SYNOPSIS:
%       [ superParticle, parameter ] = one2all(Particles, iter)
%
%   Input: 
%       Particles: Cell array of particles of size 1xN
%       iter: the number of iterations
%       outdir: output directory to store final super-particle
%
%   Output:
%       superParticle: the resulting fused particle
%       MT: Total transformation parameters (rotation+translation). MT is
%       an 4x4xNxiter matrix.
%
%   NOTE:
%       First, the function concatenates all the particles as they are.
%       Then, each particle is extracted from the stack and registered to
%       the rest. This is done until all particles are registered to the
%       rest. Once done, the whole process is iterated iter times.

% (C) Copyright 2017                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Hamidreza Heydarian, Oct 2017.

function [ superParticle, MT] = one2all(Particles, iter, oldM, outdir)

    disp('Bootstapping is started  !');
    initParticle.points = [];
    initParticle.sigma = [];
    N = numel(Particles);
   
    
    for i=1:N
        
        initParticle.points = [initParticle.points;Particles{1,i}.points];
        initParticle.sigma = [initParticle.sigma;Particles{1,i}.sigma];
        
    end
    
    superParticle{1} = initParticle.points;

    % one-to-all registration, excludes each particle from the superparticle
    % and then register it to the rest
    for j=1:iter

        tmpParticle.points = [];
        tmpParticle.sigma = [];        
        tic;
        for i=1:N

            %i = N-k+1;
            if (~mod(i,5))
                disp(['iter #' num2str(j) ' of ' num2str(iter) ': registering particle #' num2str(i) ' on initial ']);
            end
            M = Particles{1,i};
            S = delParticle(Particles, initParticle, i);
            [parameter{j,i}, ~, ~, ~, ~] = pairFitting_parallel(M, S);
            tmpParticle.points = [tmpParticle.points; transform_by_rigid2d(M.points, parameter{j,i})];
            tmpParticle.sigma = [tmpParticle.sigma; M.sigma];

        end

        a = toc;
        disp(['iter #' num2str(j) '... done in ' num2str(a) ' seconds']); 
        superParticle{j+1} = tmpParticle.points; 
        initParticle = tmpParticle;
    
    end
 
    % concatenate all previous registration parameters
    M = zeros(4,4,N,iter);
    MT = zeros(4,4,N,iter);
    for j=1:iter
        for i=1:N
            
            estAngleR = parameter{j,i}(3);
            curestAngle = wrapToPi(estAngleR);
            w_est = curestAngle * [0;0;1];
            M(1:3,1:3,i,j) = expm([0,-w_est(3),w_est(2); w_est(3),0,-w_est(1); -w_est(2),w_est(1),0]);
            M(:,4,i,j) = [parameter{j,i}(1:2)'; 0; 1]; 
%             M(:,:,i,j) = M(:,:,i,j) * oldM(:,:,i);
            R_tmp =  oldM(:,:,i);
            R_tmp(1:3,4) = 0;
            Tr_tmp = [1 0 0 -oldM(1,4,i);
                      0 1 0 -oldM(2,4,i);
                      0 0 1 -oldM(3,4,i);
                      0 0 0 oldM(4,4,i)];
            MT(:,:,i,j) = Tr_tmp * R_tmp' * M(:,:,i,j);
        end
    end

    % save to disk
    save([outdir '/superParticle'], 'superParticle');
    
    disp('Bootstapping is done  !');
    
end
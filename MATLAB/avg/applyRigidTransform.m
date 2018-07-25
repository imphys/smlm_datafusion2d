%applyRigidTransform   Applies array of rigid transform matrix to an array
% of particles
%
% SYNOPSIS:
%  newParticles = applyRigidTransform(oldParticles, M)
%
% NOTE:
%  oldParticles is an array of k cells. Each cell is a structure with fields:
%  points ans sigma. oldParticles{1,i}.points is an Nx2 matrix containing x
%  y coordinates. oldParticles{1,i}.sigma is a Nx1 matrix containing
%  localizations uncertainties. 
%  M is a 4x4xk matrix. M(1:3,1:3,i) is the rotation submatrix, M(1:2,4,i)
%  is the translation vector.

% (C) Copyright 2017                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Hamidreza Heydarian, Oct 2017.

function newParticles = applyRigidTransform(oldParticles, M)

    nParticles = numel(oldParticles);
    newParticles = cell(1,nParticles);
    flag = isfield(oldParticles{1,1}, 'points');   % 1/0 for experiment/simulation
    
    if flag
        
        % experimental data
        for i = 1:nParticles

             par = oldParticles{1,i}.points;
             tmpParticle = (par -repmat(M(1:2,4,i)',size(par,1),1))* M(1:2,1:2,i)';
             newParticles{1,i} = oldParticles{1,i};
             newParticles{1,i}.points = tmpParticle;

        end
        
    else
        
        % simulated data
        for i = 1:nParticles

             par = oldParticles{1,i}.coords(:,1:2);
             tmpParticle = (par -repmat(M(1:2,4,i)',size(par,1),1))* M(1:2,1:2,i)';
             newParticles{1,i} = oldParticles{1,i};
             newParticles{1,i}.coords(:,1:2) = tmpParticle;

        end
        
    end
    
end

% % sanitiy check
% newParticles = particles;
% nParticles = numel(newParticles);
% sr = [];
% for i=1:nParticles
%     sr = [sr;newParticles{1,i}.points];
% end
% visualizeCloud2D(sr,600,0.6,0,1);
% nParticles = numel(newParticles);
% sr = [];
% for i=1:nParticles
%     sr = [sr;newParticles{1,i}.coords(:,1:2)];
% end
% visualizeCloud2D(sr,600,0.6,0,1);
%one2all   remove localization data of index idx from fused particles fusedParticles
% of particles
%
%   SYNOPSIS:
%       S = delParticle(Particles, fusedParticles, idx)
%
%   Input: 
%       Particles: Cell array of original particles of size 1xN
%       fusedParticles: Cell array of fused particles
%       idx: index of the particle to be removed
%
%   Output:
%       S: Cell array of fusedParticles without Particles{1,idx}
%
%   NOTE:
%       This function removes the Particles{1,idx} from the stacked
%       fusedParticles.

% (C) Copyright 2017                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Hamidreza Heydarian, Oct 2017.

function S = delParticle(Particles, fusedParticles, idx)
    
    S = fusedParticles;
    particlesSize = 0;
    
    for i=1:idx
            
        particlesSize = particlesSize + numel(Particles{1,i}.points);
        
    end
    
   curParticleSize = numel(Particles{1,idx}.points);
   
   S.points(particlesSize+1:curParticleSize,:) = [];
   S.sigma(particlesSize+1:curParticleSize) = [];
   
end
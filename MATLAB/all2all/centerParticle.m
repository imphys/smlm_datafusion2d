% centerParticle   center the input particle
%
% SYNOPSIS:
%   newParticle = centerParticle(particle)
%
% INPUT
%   particle
%       The input particle
%
% OUTPUT
%   newParticle
%       The output particle
%
% Author: Hamidreza Heydarian, 2017

function newParticle = centerParticle(particle)

    cnt = mean(particle);
    newParticle = particle - repmat(cnt,size(particle,1),1);

end
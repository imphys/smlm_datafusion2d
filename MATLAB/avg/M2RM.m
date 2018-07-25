% M2RM   computes the relative transformation matrix RM from the absolute M
%
% SYNOPSIS:
%   RM = M2RM(M)
%
% INPUT
%   M
%       Absolute motion matrix (rotation+translation)
%
% OUTPUT
%   RM
%       Relative motion matrix (rotation+translation)
%
% NOTES
%   M_ij = M_j^-1 x M_i
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Hamidreza Heydarian, 2017

function RM = M2RM(M)
    
    RM = zeros(4,4,5);
    nParticles = size(M,3);
    kk = 1;
    for i=1:nParticles-1
        for j=i+1:nParticles

            relAngleR(kk) = atan2(-M(2,1,i),M(1,1,i)) + atan2(M(2,1,j),M(1,1,j));
            curestAngle = wrapToPi(relAngleR(kk));
            w_est = curestAngle * [0;0;1];
            RM(1:3,1:3,kk) = expm([0,-w_est(3),w_est(2); w_est(3),0,-w_est(1); -w_est(2),w_est(1),0]); 
            translation = -M(1:2,4,i) + M(1:2,4,j);
            RM(:,4,kk) = [translation; 0; 1]; 
            kk = kk+1;

        end
    end
            
end
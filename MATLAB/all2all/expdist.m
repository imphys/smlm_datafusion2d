% expdist           Apply rigid transform to the point set and invoke 
%                   executive to evaluate Bhattacharya cost function
%
% SYNOPSIS:
%   D = expdist(S, M, angle, t1, t2)
%
% PARAMETERS:
%   S
%      The first particle (struct with points and sigma as fields)
%   M
%      The second particle (struct with points and sigma as fields)
%   angle
%      in-plane rotation angle to be applied to M
%   [t1, t2]
%      2D translation paramteres 
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Hamidreza Heydarian, Feb 2017

function D = expdist(S, M, angle, t1, t2)

    R = [cos(angle) -sin(angle); sin(angle) cos(angle)];    % rotation matrix
    t = [t1 t2]; % translation vector
    
    % transform the model
    Mt.points = M.points * R' + repmat(t, size(M.points,1), 1);
    Mt.sigma = M.sigma;
    
    % compute the Bhatacharya cost function
    % run GPU version if it exists otherwise fall back to CPU version
    if exist('mex_expdist','file')
        D = mex_expdist(Mt.points, S.points, Mt.sigma, S.sigma);
    elseif exist('mex_expdist_cpu','file')
        D = mex_expdist_cpu(Mt.points, S.points, Mt.sigma, S.sigma);
    else
        message = ['No compiled modules found for ExpDist.\n' ...
            'Please run make in the top-level directory first.']; 
        message_id = 'MATLAB:MEXNotFound';
        error (message_id, message);


end



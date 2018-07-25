% GaussTransform   computes the Guass transform between point set A and B
% with the scale parameter scale
%
% SYNOPSIS:
%   [f, grad] = GaussTransform(A, B, scale)
%
% INPUT
%   A
%       The first point set
%
%   B
%       The second point set
%
% OUTPUT
%   f 
%       The gausstransform output. 
%
%   g
%       The gradient vector
%
% NOTES
%       The inner product between two spherical Gaussian mixtures computed 
%       using the Gauss Transform. The centers of the two mixtures are 
%       given in terms of two point sets A and B (of same dimension d)
%       represented by an mxd matrix and an nxd matrix, respectively.
%       It is assumed that all the components have the same covariance 
%       matrix represented by a scale parameter (scale).  Also, in each 
%       mixture, all the components are equally weighted.
%
% Author: bing.jian 
% Date: 2009-02-10 02:13:49 -0500 (Tue, 10 Feb 2009) 
% Revision: 121 
% Modified: Hamidreza Heydarian, 2017

function [f,g] = GaussTransform(A, B, scale)	

    % run GPU version if it exists otherwise fall back to CPU version
    if exist('mex_gausstransform','file')
        [f,g] = mex_gausstransform(A',B',scale);
        g = g';
    elseif exist('mex_gausstransform_cpu','file')
        [f,g] = mex_gausstransform_cpu(A',B',scale);
        g = g';
    else
        message = ['No compiled modules found for GaussTransform.\n' ...
            'Please run make in the top-level directory first.'];
        message_id = 'MATLAB:MEXNotFound';
        error (message_id, message);

    
end

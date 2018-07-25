% transform_by_rigid2d   perform a 2D rigid transform on a pointset and
% return the transformed pointset
%
% SYNOPSIS:
%   [result] = transform_by_rigid2d(pointset, param)
%
% INPUT
%   pointset
%       The point set to be transformed
%
%   param
%       Transformation parameters
%
% OUTPUT
%   result 
%       Transformed point set 
%
% NOTES
%       Note that here 2D rigid transform is parametrized by [translation_x,
%       translation_y, rotation_angle].
%
% RCSfile: transform_by_rigid2d.m,v 
% Author: bing.jian 
% Date: 2008-11-30 18:09:59 -0500 (Sun, 30 Nov 2008) 
% Revision: 116 
% Modified: Hamidreza Heydarian, 2017

function [result] = transform_by_rigid2d(pointset, param)

    n = size(pointset,1);
    
    % parameters of 2D-rigid transform
    delta_x = param(1);   % translation in x-dimension
    delta_y = param(2);   % translation in y-dimension
    d_theta = param(3);   % rotation in radius
    
    % form the rotation matrix
    r = [cos(d_theta) -sin(d_theta);
         sin(d_theta)  cos(d_theta)];
     
    % first rotate
    result(:,1:2) = pointset(:,1:2) * r' ;
    
    % then translate
    result(:,1:2) = result(:,1:2) + ones(n,1)*[delta_x delta_y];
    
end
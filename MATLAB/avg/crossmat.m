% crossmat   Cross-product matrix (anti-symmetric)
%
% SYNOPSIS:
%   M=crossmat(x)
%
% INPUT
%   x 
%       input vector 1x3
%
% OUTPUT
%   M
%       cross product matrix
%
% NOTES
%       Given vector x, finds M such that cross(x,y)=My
%
% Author: Avishek Chatterjee, Venu Madhav Govindu
% 

function M=crossmat(x)

    M=[0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];
    
end
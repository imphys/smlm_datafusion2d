% R2w   Rotation matrix to vector
%
% SYNOPSIS:
%   [w]=R2w(R)
%
% INPUT
%   R 
%       3x3 rotation matrix
%
% OUTPUT
%   w
%       3 dimensional matrix such that R=expm(crossmat(w))
%
% NOTES
%       Kanatani : Geometric Computation for Machine Vision
%       Page 102 Eqn 5.7. A 3D rotation can be represented as a 3 element 
%       vector ! such that the angle of the rotation  is equal to the 
%       magnitude of the vector, and the axis of rotation n is parallel to
%       the direction of the vector [16].
%
% Author: Avishek Chatterjee, Venu Madhav Govindu
% 

function [w]=R2w(R)

    w=[R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]/2;
    s=norm(w);
    if(s)
        w=w/s*atan2(s,(trace(R)-1)/2);
    end
    
end


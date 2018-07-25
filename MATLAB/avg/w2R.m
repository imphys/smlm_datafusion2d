% w2R   Vector to rotation matrix
%
% SYNOPSIS:
%   R=w2R(w)
%
% INPUT
%   w
%       3 dimensional matrix such that R=expm(crossmat(w))
%
% OUTPUT
%   R 
%       3x3 rotation matrix
%
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

function R=w2R(w)

    omega=norm(w);
    if(omega)
        n=w/omega;

        s=sin(omega);
        c=cos(omega);
        cc=1-c;

        n1=n(1);                n2=n(2);                n3=n(3);
        n12cc=n1*n2*cc;         n23cc=n2*n3*cc;         n31cc=n3*n1*cc;
        n1s=n1*s;               n2s=n2*s;               n3s=n3*s;

        R(1,1)=c+n1*n1*cc;      R(1,2)=n12cc-n3s;       R(1,3)=n31cc+n2s;
        R(2,1)=n12cc+n3s;       R(2,2)=c+n2*n2*cc;      R(2,3)=n23cc-n1s;
        R(3,1)=n31cc-n2s;       R(3,2)=n23cc+n1s;       R(3,3)=c+n3*n3*cc;
    else
        R=eye(3);
    end
    
end
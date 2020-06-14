%connectGraph   This function keeps increasing the threshold on the error
%to include more connections, untill the graph is connected.
%
%   SYNOPSIS:
%       [error_idx] = connectGraph(I,error)
%
%   Input: 
%       I: 2xN matrix of indices indications the connections
%       error: vector of error values for every combinations in I
%
%   Output:
%       error_idx: indices of the to-be-removed connection in I
%
%   NOTE:
%
% (C) Copyright 2017                    QI Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Teun Huijben, June 2020.

function error_idx = connectGraph(I,error)

    [a,b] = sort(abs(error));           %sort the errors
    first_above5 = find(a>5,1);         %find the index of lowest error

    error_idx = abs(error)>5;           %indices of all errors
    I_new = I; 
    I_new(:,error_idx) = [];            %remove errors from graph

    start = first_above5; 
    while (isConnected(I_new)==0 | max(I_new(:))<max(I(:)))
        error_idx = (abs(error)>abs(error(b(start))));  %indices of all errors including new one
        I_new = I; 
        I_new(:,error_idx) = []; 
        start = start+1; 
    end
 
end

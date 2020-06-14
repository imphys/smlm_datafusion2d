%isConnected   This function checks whether the graph is connected
%
%   SYNOPSIS:
%       [connected] = isConnected(I)
%
%   Input: 
%       I: 2xN matrix of indices indications the connections
%
%   Output:
%       connected: Boolean indicating whether is connected
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

function connected = isConnected(I)

    Gr = graph(I(1,:),I(2,:));
    clusters = conncomp(Gr);
    if max(clusters)==1
        connected = true;
    else
        connected = false; 
    end

end

